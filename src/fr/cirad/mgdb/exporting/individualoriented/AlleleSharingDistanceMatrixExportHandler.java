package fr.cirad.mgdb.exporting.individualoriented;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.EigenstratExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.dist.AlleleSharingDistanceMatrixCalculator;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ExperimentalFeature;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class AlleleSharingDistanceMatrixExportHandler.
 * Computes ASD (Allele Sharing Distance) matrix from Eigenstrat genotype data.
 */
public class AlleleSharingDistanceMatrixExportHandler extends EigenstratExportHandler implements ExperimentalFeature {

    private int nMaxMissingDataPercentageForIndividuals = 50;
    
    private final static Pattern missingAllelePattern = Pattern.compile("[9]"); // Eigenstrat uses 9 for missing
    
    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(AlleleSharingDistanceMatrixExportHandler.class);

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
    }

    public AlleleSharingDistanceMatrixExportHandler() {
    }
    
    public AlleleSharingDistanceMatrixExportHandler(int nMaxMissingDataPercentageForIndividuals) {
        super();
        this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }
    
    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "ASD-DIST";
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports a zipped distance matrix file calculated according to Allele Sharing Distance (ASD) based on an Eigenstrat format export. Individuals with more than " + nMaxMissingDataPercentageForIndividuals + "% missing data are automatically excluded.";
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getSupportedVariantTypes()
     */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }
    
    @Override
    public String getExportArchiveExtension() {
        return "zip";
    }

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        
        // Build sorted individual list
        TreeSet<String> indSet = new TreeSet<>(new AlphaNumericComparator<String>());
        for (Callset cs : callSetsToExport)
            indSet.add(workWithSamples ? cs.getSampleId() : cs.getIndividual());
        List<String> sortedIndividuals = new ArrayList<>(indSet);

        // Build individual positions map
        Map<String, Integer> individualPositions = IExportHandler.buildIndividualPositions(callSetsToExport, workWithSamples);

        if (progress.isAborted())
            return;

        // Build variant query
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : 
            (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : 
                (varQueryWrapper.getBareQueries().iterator().hasNext() ? 
                    new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

        // Save warnings to temp file
        File warningFile = File.createTempFile("export_warnings_", "");
        File tempGenoFile = File.createTempFile("eigenstrat_for_asd__", "");
        
        try {
            // Write Eigenstrat genotype file
            OutputStream genoOS = new FileOutputStream(tempGenoFile);
            ExportOutputs exportOutputs = writeGenotypeFile(genoOS, sModule, nAssemblyId, workWithSamples, 
                annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, 
                markerCount, markerSynonyms, callSetsToExport, individualsByPop);
            genoOS.close();

            // Collect warnings from export
            FileOutputStream warningOS = new FileOutputStream(warningFile);
            for (File f : exportOutputs.getWarningFiles()) {
                if (f != null && f.length() > 0) {
                    BufferedReader in = new BufferedReader(new FileReader(f));
                    String sLine;
                    while ((sLine = in.readLine()) != null)
                        warningOS.write((sLine + "\n").getBytes());
                    in.close();
                    f.delete();
                }
            }
            warningOS.close();

            progress.moveToNextStep();

            // Process genotype data from temp file
            List<String> filteredIndividuals = new ArrayList<>();
            StringBuilder missingDataWarnings = new StringBuilder();
            
            // Read and filter individuals based on missing data
            try (BufferedReader reader = new BufferedReader(new FileReader(tempGenoFile))) {
                String line;
                int individualIndex = 0;
                
                while ((line = reader.readLine()) != null && individualIndex < sortedIndividuals.size()) {
                    String individualName = sortedIndividuals.get(individualIndex);
                    String genotypeLine = line.trim();
                    
                    // Calculate missing data percentage (9 = missing in Eigenstrat)
                    Matcher matcher = missingAllelePattern.matcher(genotypeLine);
                    int missingAlleleCount = 0;
                    while (matcher.find())
                        missingAlleleCount++;
                    
                    int missingPercent = genotypeLine.length() > 0 ? 
                        missingAlleleCount * 100 / genotypeLine.length() : 100;
                    
                    if (missingPercent > nMaxMissingDataPercentageForIndividuals)
                        missingDataWarnings.append("- Excluding individual ").append(individualName).append(" from ASD export, it has too much missing data: ").append(missingPercent).append("%\n");
                    else
                        filteredIndividuals.add(individualName);
                    
                    individualIndex++;
                    
                    // Update progress
                    if (individualIndex % 100 == 0 && progress != null) {
                        int pct = (int) ((individualIndex / (double) sortedIndividuals.size()) * 50); // First half of progress
                        progress.setCurrentStepProgress(pct);
                        
                        if (progress.getError() != null || progress.isAborted()) {
                            LOG.debug("ASD export aborted during genotype processing");
                            return;
                        }
                    }
                }
            }

            if (filteredIndividuals.size() < 2) {
                throw new Exception("Not enough individuals with sufficient data for ASD calculation. Only " + 
                    filteredIndividuals.size() + " individuals passed the " + nMaxMissingDataPercentageForIndividuals + "% missing data threshold.");
            }

            LOG.debug("Proceeding with " + filteredIndividuals.size() + " individuals out of " + 
                sortedIndividuals.size() + " for ASD calculation");

            // Reread genotype data for filtered individuals
            progress.setCurrentStepProgress(50);
            
            // Convert to byte arrays for faster processing
            final int numMarkers = (int) markerCount;
            byte[][] genotypeBytes = readFilteredGenotypeData(tempGenoFile, sortedIndividuals, filteredIndividuals, numMarkers, progress);

            if (genotypeBytes == null || progress.isAborted() || progress.getError() != null) {
                return;
            }

            progress.moveToNextStep();

            // Compute ASD distance matrix
            double[][] distanceMatrix = new AlleleSharingDistanceMatrixCalculator(genotypeBytes/*, filteredIndividuals.size()*/, numMarkers).calculate(progress);
            

            if (distanceMatrix == null || progress.isAborted() || progress.getError() != null) {
                return;
            }

            // Create export name
            Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
            String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + 
                "__" + markerCount + "variants__" + filteredIndividuals.size() + "individuals";

            // Create ZIP output
            ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
            
            // Add metadata if requested
            if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty()) {
                IExportHandler.addMetadataEntryIfAny(sModule + "__" + filteredIndividuals.size() + 
                    (workWithSamples ? "sample" : "individual") + "s_metadata.tsv", sModule, sExportingUser, 
                    new TreeSet<>(filteredIndividuals), individualMetadataFieldsToExport, zos, 
                    (workWithSamples ? "sample" : "individual"), workWithSamples);
            }

            // Write distance matrix
            finalizeExportUsingDistanceMatrix(filteredIndividuals, exportName, distanceMatrix, zos, progress);

            // Add warnings if any
            if (progress.getError() == null && !progress.isAborted()) {
                if (warningFile.length() > 0 || missingDataWarnings.length() > 0) {
                    progress.addStep("Adding lines to warning file");
                    progress.moveToNextStep();
                    progress.setPercentageEnabled(false);
                    zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));

                    if (missingDataWarnings.length() > 0)
                        zos.write(missingDataWarnings.toString().getBytes());

                    if (warningFile.length() > 0) {
                        try (BufferedReader in = new BufferedReader(new FileReader(warningFile))) {
                            String sLine;
                            int nWarningCount = 0;
                            while ((sLine = in.readLine()) != null) {
                                zos.write((sLine + "\n").getBytes());
                                progress.setCurrentStepProgress(nWarningCount++);
                            }
                            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
                        }
                    }

                    zos.closeEntry();
                }
            }

            zos.finish();
            zos.close();
        }
        catch (OutOfMemoryError oome) {
            LOG.error("Not enough RAM for ASD calculation", oome);
            progress.setError("Unable to compute ASD (dataset is too large): " + oome.getMessage());
        }
        finally {
            warningFile.delete();
            tempGenoFile.delete();
        }

        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }
    
    /**
     * Reads genotype data from a file into byte arrays, producing a matrix of shape [nIndividuals][nMarkers].
     * Supports filtering individuals by name.
     *
     * File format:
     * - Each row = a SNP / marker
     * - Each column = an individual
     * - Characters = '0', '1', '2', '9' (missing)
     *
     * @param genoFile          input file
     * @param allIndividuals    list of all individual names in file order
     * @param filteredIndividuals subset to keep (may be same as allIndividuals)
     * @param numMarkers        number of markers / rows in the file
     * @param progress          progress indicator
     * @return [nFilteredIndividuals][numMarkers] numeric byte matrix
     */
    public static byte[][] readFilteredGenotypeData(
            File genoFile,
            List<String> allIndividuals,
            List<String> filteredIndividuals,
            int numMarkers,
            ProgressIndicator progress) throws IOException {

        int nAll = allIndividuals.size();
        int nFiltered = filteredIndividuals.size();

        // Map individual name → column index
        Map<String, Integer> individualIndexMap = new HashMap<>();
        for (int i = 0; i < nAll; i++) {
            individualIndexMap.put(allIndividuals.get(i), i);
        }

        // Prepare filtered indices (columns)
        int[] filteredIndices = new int[nFiltered];
        for (int i = 0; i < nFiltered; i++) {
            filteredIndices[i] = individualIndexMap.get(filteredIndividuals.get(i));
        }

        // Allocate output: [nFilteredIndividuals][numMarkers]
        byte[][] genotypeBytes = new byte[nFiltered][numMarkers];

        try (BufferedReader reader = new BufferedReader(new FileReader(genoFile))) {
            String line;
            int markerIndex = 0;

            while ((line = reader.readLine()) != null && markerIndex < numMarkers) {
                line = line.trim();
                if (line.length() < nAll) {
                    throw new IllegalStateException("Line " + markerIndex + " too short: expected " + nAll + " chars");
                }

                // Copy only filtered individuals
                for (int fi = 0; fi < nFiltered; fi++) {
                    int col = filteredIndices[fi];
                    char c = line.charAt(col);
                    if (c >= '0' && c <= '2') {
                        genotypeBytes[fi][markerIndex] = (byte) (c - '0');
                    } else {
                        genotypeBytes[fi][markerIndex] = 9; // missing
                    }
                }

                markerIndex++;

                // Progress update
                if (progress != null && markerIndex % 100 == 0) {
                    int pct = 50 + (int) ((markerIndex / (double) numMarkers) * 25);
                    progress.setCurrentStepProgress(pct);

                    if (progress.isAborted() || progress.getError() != null) {
                        System.out.println("ASD export aborted during data reading");
                        return null;
                    }
                }
            }

            if (markerIndex != numMarkers) {
                throw new IllegalStateException("Expected " + numMarkers + " markers, but got " + markerIndex);
            }
        }

        return genotypeBytes;
    }

    
//    /**
//     * Calculate ASD distance matrix using parallel processing.
//     */
//    private double[][] calculateDistanceMatrix(byte[][] genotypeBytes, int nIndividuals, int nMarkers, ProgressIndicator progress) throws InterruptedException {
//        
//        // Upper-triangle allocation (saves 50% memory)
//        final double[][] dist = new double[nIndividuals][];
//        for (int i = 0; i < nIndividuals; i++)
//            dist[i] = new double[nIndividuals - i - 1];
//
//        final long totalPairs = (long) nIndividuals * (nIndividuals - 1) / 2;
//        final LongAdder completed = new LongAdder();
//
//        final int cpu = Runtime.getRuntime().availableProcessors();
//        final int threads = Math.max(1, cpu / 2);
//        
//        LOG.debug("Launching ASD calculation on " + threads + " threads for " + nIndividuals + " individuals with " + nMarkers + " markers (" + totalPairs + " pairs)");
//
//        ExecutorService pool = Executors.newFixedThreadPool(threads);
//        List<Future<Void>> futures = new ArrayList<>();
//
//        // Process in blocks of individuals
//        final int BLOCK_SIZE = 32;
//        List<int[]> blocks = new ArrayList<>();
//        for (int start = 0; start < nIndividuals; start += BLOCK_SIZE) {
//            int from = start;
//            int to = Math.min(nIndividuals, start + BLOCK_SIZE);
//            blocks.add(new int[]{from, to});
//        }
//
//        for (int[] block : blocks) {
//            final int from = block[0];
//            final int to = block[1];
//
//            futures.add(pool.submit(() -> {
//                long localCompleted = 0;
//                final long PROGRESS_STEP = 50_000;
//
//                for (int i = from; i < to; i++) {
//                    // Check for cancellation only once per row
//                    if (progress != null && (progress.getError() != null || progress.isAborted()))
//                        return null;
//                    
//                    final byte[] genoI = genotypeBytes[i];
//
//                    for (int j = i + 1; j < nIndividuals; j++) {
//                        double asd = calculateASD(genoI, genotypeBytes[j], nMarkers);
//                        dist[i][j - i - 1] = asd;
//                        localCompleted++;
//                    }
//                    
//                    // Update shared counter after completing each row to reduce contention
//                    completed.add(localCompleted);
//                    long done = completed.sum();
//                    
//                    // Throttle progress updates
//                    if (done % PROGRESS_STEP < (nIndividuals - i - 1) && progress != null) {
//                        int pct = 75 + (int) ((done / (double) totalPairs) * 25); // Last 25% of progress
//                        progress.setCurrentStepProgress(pct);
//                    }
//                    
//                    localCompleted = 0;
//                }
//                
//                // Add any remaining local count
//                if (localCompleted > 0) {
//                    completed.add(localCompleted);
//                }
//                
//                return null;
//            }));
//        }
//
//        pool.shutdown();
//        
//        // Wait with regular cancellation checks
//        while (!pool.awaitTermination(1, TimeUnit.SECONDS)) {
//            if (progress != null && (progress.getError() != null || progress.isAborted())) {
//                pool.shutdownNow();
//                return null;
//            }
//        }
//
//        if (progress != null && (progress.getError() != null || progress.isAborted()))
//            return null;
//
//        // Wait for all futures to complete
//        for (Future<Void> f : futures) {
//            try {
//                f.get();
//            } catch (ExecutionException e) {
//                if (progress != null)
//                    progress.setError("ASD computation failed: " + e.getCause().getMessage());
//                throw new RuntimeException("ASD computation failed", e.getCause());
//            }
//        }
//
//        if (progress != null)
//            progress.setCurrentStepProgress(100);
//        
//        LOG.debug("ASD matrix calculation completed");
//        return dist;
//    }
    
	protected void finalizeExportUsingDistanceMatrix(List<String> sequenceNames, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
	    zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
	    java.text.NumberFormat formatter = new java.text.DecimalFormat("#0.000000", new DecimalFormatSymbols(Locale.US)); 
	    
	    int n = sequenceNames.size();	    
	    zos.write(("\t" + n).getBytes());	// Write header with number of individuals
	    
	    // Write complete distance matrix (not just upper triangle)
	    for (int i = 0; i < n; i++) {
	        StringBuilder rowBuilder = new StringBuilder();
	        rowBuilder.append("\n").append(sequenceNames.get(i));
	        
	        for (int j = 0; j < n; j++) {
	            double distance;
	            if (i == j)
	                distance = 0.0;  // Diagonal: distance to self is 0
	            else if (i < j)
	            	distance = distanceMatrix[i][j - i - 1];	// Upper triangle: retrieve from distanceMatrix[i][j-i-1]
	            else
	                distance = distanceMatrix[j][i - j - 1];	// Lower triangle: symmetric to upper triangle
  
	            rowBuilder.append(" ").append(formatter.format(distance));
	        }
	        
	        zos.write(rowBuilder.toString().getBytes());
	        
	        if (n > 1000 && i % 500 == 0 && progress != null)	// Update progress for large datasets
	            progress.setCurrentStepProgress((int) ((i / (double) n) * 100));
	    }
	    
	    zos.closeEntry();
	}


    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{
            "Exporting to Eigenstrat format", 
            "Filtering individuals with missing data",
            "Calculating Allele Sharing Distance matrix"
        });
    }

    @Override
    public String[] getExportDataFileExtensions() {
        return new String[] {"mtx"};
    }

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
    
//	 // Lookup table for shared alleles between genotypes
//	 // Index calculation: genotype1 * 10 + genotype2
//	 // Values 0-2 represent hom ref, het, hom alt; 9 represents missing
//	 private static final byte[] SHARED_ALLELES_LOOKUP = new byte[100];
//	
//	 static {
//	     // Initialize all to -1 (invalid)
//	     Arrays.fill(SHARED_ALLELES_LOOKUP, (byte) -1);
//	     
//	     // Valid comparisons (symmetric, so we fill both directions)
//	     // 0 vs 0: 2 shared
//	     SHARED_ALLELES_LOOKUP[0 * 10 + 0] = 2;
//	     
//	     // 0 vs 1: 1 shared
//	     SHARED_ALLELES_LOOKUP[0 * 10 + 1] = 1;
//	     SHARED_ALLELES_LOOKUP[1 * 10 + 0] = 1;
//	     
//	     // 0 vs 2: 0 shared
//	     SHARED_ALLELES_LOOKUP[0 * 10 + 2] = 0;
//	     SHARED_ALLELES_LOOKUP[2 * 10 + 0] = 0;
//	     
//	     // 1 vs 1: 2 shared
//	     SHARED_ALLELES_LOOKUP[1 * 10 + 1] = 2;
//	     
//	     // 1 vs 2: 1 shared
//	     SHARED_ALLELES_LOOKUP[1 * 10 + 2] = 1;
//	     SHARED_ALLELES_LOOKUP[2 * 10 + 1] = 1;
//	     
//	     // 2 vs 2: 2 shared
//	     SHARED_ALLELES_LOOKUP[2 * 10 + 2] = 2;
//	     
//	     // 9 (missing) comparisons are left as -1 to indicate skip
//	 }
//	
//	 /**
//	  * Calculate ASD between two individuals using optimized lookup table.
//	  * 
//	  * This version:
//	  * - Uses numeric comparisons (faster than char comparisons)
//	  * - Employs lookup table for O(1) shared allele determination
//	  * - Corrects the heterozygote logic
//	  * - Eliminates redundant conditionals
//	  */
//	 private double calculateASD(byte[] geno1, byte[] geno2, int nMarkers) {
//	     int totalComparisons = 0;
//	     int sharedAlleles = 0;
//	     
//	     for (int m = 0; m < nMarkers; m++) {
//	         // Convert ASCII to numeric (assuming input is ASCII '0', '1', '2', '9')
//	         int g1 = geno1[m] - '0';  // '0' -> 0, '1' -> 1, '2' -> 2, '9' -> 9
//	         int g2 = geno2[m] - '0';
//	         
//	         // Lookup shared alleles (-1 means skip due to missing data)
//	         int lookupIndex = g1 * 10 + g2;
//	         if (lookupIndex >= 0 && lookupIndex < 100) {
//	             byte shared = SHARED_ALLELES_LOOKUP[lookupIndex];
//	             
//	             if (shared >= 0) {  // Valid comparison (not missing)
//	                 totalComparisons++;
//	                 sharedAlleles += shared;
//	             }
//	         }
//	     }
//	     
//	     if (totalComparisons == 0)
//	         return 1.0; // No comparable data → maximum distance
//	     
//	     // ASD = 1 - (shared_alleles / (2 * total_comparisons))
//	     return 1.0 - (sharedAlleles / (2.0 * totalComparisons));
//	 }
	
//	 /**
//	  * Alternative version with explicit numeric handling if bytes are already numeric.
//	  * Use this if your bytes are stored as actual numbers 0, 1, 2, 9 (not ASCII).
//	  */
//	 private double calculateASDNumeric(byte[] geno1, byte[] geno2, int nMarkers) {
//	     int totalComparisons = 0;
//	     int sharedAlleles = 0;
//	     
//	     for (int m = 0; m < nMarkers; m++) {
//	         int g1 = geno1[m];
//	         int g2 = geno2[m];
//	         
//	         // Fast missing data check
//	         if (g1 == 9 || g2 == 9)
//	             continue;
//	         
//	         totalComparisons++;
//	         
//	         // Lookup shared alleles
//	         int lookupIndex = g1 * 10 + g2;
//	         sharedAlleles += SHARED_ALLELES_LOOKUP[lookupIndex];
//	     }
//	     
//	     if (totalComparisons == 0)
//	         return 1.0;
//	     
//	     return 1.0 - (sharedAlleles / (2.0 * totalComparisons));
//	 }
//	
//	 /**
//	  * Ultra-optimized version with manual inlining.
//	  * Use this for maximum performance if profiling shows the lookup is a bottleneck.
//	  */
//	 private double calculateASDInlined(byte[] geno1, byte[] geno2, int nMarkers) {
//	     int totalComparisons = 0;
//	     int sharedAlleles = 0;
//	     
//	     for (int m = 0; m < nMarkers; m++) {
//	         int g1 = geno1[m] - '0';
//	         int g2 = geno2[m] - '0';
//	         
//	         // Skip missing data (9)
//	         if (g1 == 9 || g2 == 9)
//	             continue;
//	         
//	         totalComparisons++;
//	         
//	         // Inline the lookup logic for maximum speed
//	         // Same genotype
//	         if (g1 == g2) {
//	             sharedAlleles += 2;  // 0==0, 1==1, or 2==2 all share 2 alleles
//	         }
//	         // Different genotypes
//	         else if ((g1 == 0 && g2 == 2) || (g1 == 2 && g2 == 0)) {
//	             // Different homozygotes: 0 shared
//	             // sharedAlleles += 0; (no-op, but explicit for clarity)
//	         }
//	         else {
//	             // One or both heterozygous, different genotypes: 1 shared
//	             // This covers: 0-1, 1-0, 1-2, 2-1
//	             sharedAlleles += 1;
//	         }
//	     }
//	     
//	     if (totalComparisons == 0)
//	         return 1.0;
//	     
//	     return 1.0 - (sharedAlleles / (2.0 * totalComparisons));
//	 }
}