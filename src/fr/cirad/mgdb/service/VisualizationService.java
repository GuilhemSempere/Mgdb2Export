/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 - 2025, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.service;

import java.io.ByteArrayOutputStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.stereotype.Component;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.MongoCommandException;
import com.mongodb.client.AggregateIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.model.Projections;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.HapMapExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.MgdbChartRequest;
import fr.cirad.model.MgdbVcfFieldPlotRequest;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryBuilder;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import fr.cirad.tools.query.GroupedExecutor;
import fr.cirad.tools.query.GroupedExecutor.TaskWrapper;
import fr.cirad.tools.security.base.AbstractTokenManager;

/**
 * A service class responsible for generating dynamic visualization data
 *
 * @author sempere
 */

@Component
public class VisualizationService {
    protected static final Logger LOG = Logger.getLogger(VisualizationService.class);

    final private int MIN_INTERVAL_SIZE = 10; // minimal interval size in bp
	private int VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK = 5;

	/**
	 * Set to <= 0 to disable prechecking
	 * @param threshold
	 */
	public void setVariantsPerIntervalThresholdForPrecheck(int threshold) {
		this.VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK = threshold;
	}

	public int getVariantsPerIntervalThresholdForPrecheck() {
		return this.VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK;
	}

    protected boolean findDefaultRangeMinMax(MgdbChartRequest mdr, String tmpCollName /* if null, main variant coll is used*/) throws Exception
    {
        if (mdr.getDisplayedRangeMin() != null && mdr.getDisplayedRangeMax() != null) {
            LOG.info("findDefaultRangeMinMax skipped because min-max values already set");
            return true;    // nothing to do
        }

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(mdr.getVariantSetId());
        List<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        Long[] minMaxFound = new Long[2];   // will be filled in by the method below
        boolean retVal = Helper.findDefaultRangeMinMax(info[0], projIDs, tmpCollName, mdr.getDisplayedVariantType() == null ? null : Arrays.asList(mdr.getDisplayedVariantType()), Arrays.asList(mdr.getDisplayedSequence()), mdr.getStart(), mdr.getEnd(), minMaxFound);

        mdr.setDisplayedRangeMin(minMaxFound[0]);
        mdr.setDisplayedRangeMax(minMaxFound[1]);

        return retVal;
    }

    public Map<Long, Long> selectionDensity(MgdbChartRequest gdr, String processId) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(processId))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(processId, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "variant density on sequence " + gdr.getDisplayedSequence() });
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        final ConcurrentHashMap<Long, Long> result = new ConcurrentHashMap<Long, Long>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, nTempVarCount == 0 ? null : tmpVarColl.getNamespace().getCollectionName())) {
                progress.setError("selectionDensity: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filters).");
                return result;
            }

        final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        final ProgressIndicator finalProgress = progress;

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), nTempVarCount == 0 ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), 0L);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (!finalProgress.isAborted()) {
                        long partialCount = mongoTemplate.getCollection(nTempVarCount == 0 ? mongoTemplate.getCollectionName(VariantData.class) : tmpVarColl.getNamespace().getCollectionName()).countDocuments(intervalEntry.getValue());
                        nTotalTreatedVariantCount.addAndGet((int) partialCount);
                        result.put(intervalEntry.getKey(), partialCount);
                        finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    }
                    else {
                        for (Future<Void> f : threadsToWaitFor)
                            f.cancel(true);
                        LOG.debug("Cancelled query threads for process " + finalProgress.getProcessId());
                    }
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId()); // important to be sure that all tasks in the group are executed before the queue purges it
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionDensity treated " + nTotalTreatedVariantCount.get() + " variants on " + (nTempVarCount == 0 ? "main" : "temporary") + " collection for sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<Long, Long>(result);
    }

    // TODO: Refactor this?
    private void mergeVariantQueryDBList(BasicDBObject matchStage, BasicDBList variantQueryDBList) {
        Iterator<Object> queryItems = variantQueryDBList.iterator();
        while (queryItems.hasNext()) {
            BasicDBObject queryItem = (BasicDBObject)queryItems.next();
            for (String key : queryItem.keySet()) {
                if (queryItem.get(key) instanceof BasicDBObject) {
                    BasicDBObject queryItemElement = (BasicDBObject)queryItem.get(key);
                    if (matchStage.containsKey(key)) {
                        if (matchStage.get(key) instanceof BasicDBObject) {
                            BasicDBObject matchStageElement = (BasicDBObject)matchStage.get(key);
                            for (String elementKey : queryItemElement.keySet()) {
                                if (matchStageElement.containsKey(elementKey)) {
                                    if (elementKey.equals("$lt") || elementKey.equals("$lte")) {
                                        matchStageElement.put(elementKey, Math.min(matchStageElement.getLong(elementKey), queryItemElement.getLong(elementKey)));
                                    } else if (elementKey.equals("$gt") || elementKey.equals("$gte")) {
                                        matchStageElement.put(elementKey, Math.max(matchStageElement.getLong(elementKey), queryItemElement.getLong(elementKey)));
                                    } else {
                                        matchStageElement.put(elementKey, queryItemElement.get(elementKey));
                                    }
                                } else {
                                    matchStageElement.put(elementKey, queryItemElement.get(elementKey));
                                }
                            }
                        } else {
                            matchStage.put(key, queryItemElement);
                        }
                    } else {
                        matchStage.put(key, queryItemElement);
                    }
                } else {
                    matchStage.put(key, queryItem.get(key));
                }
            }
        }
    }

    private HashMap<Long, BasicDBObject> getIntervalQueries(
            int nIntervalCount,
            Collection<String> sequences,
            String variantType,
            long rangeMin,
            long rangeMax,
            BasicDBList variantQueryDBListToMerge,
            MongoCollection<Document> precheckCollection // can be null if we don't want to precheck
    ) throws InterruptedException, ExecutionException {

        HashMap<Long, BasicDBObject> result = new HashMap<>();
        int rangeSize = (int) (rangeMax - rangeMin);
        if (rangeSize < 0)
        	return result;

        int actualIntervalCount = Math.min(rangeSize, Math.max(1, nIntervalCount));
        long intervalSize = (long) Math.ceil((double) rangeSize / actualIntervalCount);

        // if there are really few variants per interval on average we consider it is worth pre-checking for empty intervals
        String refPosPathWithDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
        ExecutorService preCheckExecutor = VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK <= 0 || precheckCollection == null || (float) precheckCollection.countDocuments(new BasicDBObject(refPosPathWithDot + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$gte", rangeMin).append("$lte", rangeMax))) / actualIntervalCount > VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK  ? null : Executors.newFixedThreadPool(Math.max(1, Runtime.getRuntime().availableProcessors() / 2));
        AtomicInteger keptQueryCount = new AtomicInteger(0);

        List<CompletableFuture<Void>> futures = new ArrayList<>();
        for (int i = 0; i < actualIntervalCount; i++) {
            long start = rangeMin + i * intervalSize;
            long end = (i == actualIntervalCount - 1) ? rangeMax : start + intervalSize;

            BasicDBObject query = new BasicDBObject();
            if (sequences != null && !sequences.isEmpty())
                query.put(refPosPathWithDot + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", sequences));
            if (variantType != null)
                query.put(VariantData.FIELDNAME_TYPE, variantType);

            BasicDBObject positionSettings = new BasicDBObject("$gte", start);
            positionSettings.put(i < actualIntervalCount - 1 ? "$lt" : "$lte", end);
            query.put(refPosPathWithDot + ReferencePosition.FIELDNAME_START_SITE, positionSettings);

            // Async precheck
            if (preCheckExecutor != null) {
                CompletableFuture<Void> future = CompletableFuture.runAsync(() -> {
                    boolean hasData = precheckCollection.find(query)
                            .projection(Projections.include("_id"))
                            .limit(1)
                            .iterator()
                            .hasNext();
                    synchronized (result) {
                        if (hasData) {
                            if (variantQueryDBListToMerge != null && !variantQueryDBListToMerge.isEmpty())
                                mergeVariantQueryDBList(query, variantQueryDBListToMerge);
                            result.put(start, query);
                            keptQueryCount.incrementAndGet();
                        } else {
                            result.put(start, null); // keep interval accounted for
                        }
                    }
                }, preCheckExecutor);
                futures.add(future);
            } else {
                if (variantQueryDBListToMerge != null && !variantQueryDBListToMerge.isEmpty())
                    mergeVariantQueryDBList(query, variantQueryDBListToMerge);
                result.put(start, query);
            }
        }

        if (!futures.isEmpty()) {
            CompletableFuture.allOf(futures.toArray(new CompletableFuture[0])).get();
            preCheckExecutor.shutdown();
            if (keptQueryCount.get() < actualIntervalCount)
            	LOG.debug("Precheck led to reduce number of interval queries from " + actualIntervalCount + " to " + keptQueryCount.get());
        }

        return result;
    }

    public Map<Long, Double> selectionFst(MgdbChartRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Fst estimate on sequence " + gdr.getDisplayedSequence()});
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);
        Collection<BasicDBList> variantRunDataQueries = varQueryWrapper.getVariantRunDataQueries();
        final BasicDBList variantQueryDBList = variantRunDataQueries.size() == 1 ? variantRunDataQueries.iterator().next() : new BasicDBList();

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        final boolean useTempColl = (nTempVarCount != 0);
        final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : mongoTemplate.getCollectionName(VariantRunData.class);
        final ConcurrentHashMap<Long, Double> result = new ConcurrentHashMap<Long, Double>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, useTempColl ? usedVarCollName : null)) {
                progress.setError("selectionFst: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filters).");
                return result;
            }

        final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);
        final ProgressIndicator finalProgress = progress;

        List<BasicDBObject> baseQuery = buildFstQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

            List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            if (!useTempColl)
	            intervalEntry.getValue().append(VariantData.FIELDNAME_KNOWN_ALLELES + ".2", new BasicDBObject("$exists", false)); // exclude multi-allelic variants from MAF calculation
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
            Thread t = new Thread() {
                public void run() {
                    if (finalProgress.isAborted())
                        return;

                    try {
                        Iterator<Document> it = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(true).iterator();

                        if (finalProgress.isAborted())
                            return;

                        /* Structure of a resulting document : {
                         *      _id: ...,
                         *      alleleMax: 1,
                         *      populations: [
                         *          {sampleSize: 100, alleles: [
                         *              {allele: 0, alleleFrequency: 0.45, heterozygoteFrequency: 0.31},
                         *              {allele: 1, alleleFrequency: 0.55, heterozygoteFrequency: 0.31},
                         *          ]},
                         *          {...}
                         *      ]
                         * }
                         */

                        double weightedFstSum = 0;
                        double fstWeight = 0;
                        int partialCount = 0;

                        while (it.hasNext()) {
                            partialCount++;
                            Document variantResult = it.next();

                            List<Document> populations = variantResult.getList(FST_RES_POPULATIONS, Document.class);
                            if (populations.size() < 2) {
                                // Can not compute Fst with a single population
                                // One of the populations has no valid data
                                continue;
                            }
                            int numPopulations = populations.size();  // r : Number of samples to consider
                            int numAlleles = variantResult.getInteger(FST_RES_ALLELEMAX) + 1;

                            // Transposition to [allele][sample] instead of the original [sample][allele] is important to simplify further computations
                            int[] sampleSizes = new int[numPopulations];  // n_i = sampleSizes[population] : Size of the population samples (with missing data filtered out)
                            double[][] alleleFrequencies = new double[numAlleles][numPopulations];  // p_i = alleleFrequencies[allele][population] : Allele frequency in the given population
                            double[][] hetFrequencies = new double[numAlleles][numPopulations];  // h_i = hetFrequencies[allele][population] : Proportion of heterozygotes with the given allele in the given population
                            //double[] averageAlleleFrequencies = new double[numAlleles];  // p¯ = averageAlleleFrequencies[allele] : Average frequency of the allele over all populations
                            //double[] alleleVariance = new double[numAlleles];  // s² = alleleVariance[allele] : Variance of the allele frequency over the populations
                            //double[] averageHetFrequencies = new double[numAlleles];  // h¯ = averageHetFrequencies : Proportion of heterozygotes with the given allele over all populations

                            //Arrays.fill(averageAlleleFrequencies, 0);
                            //Arrays.fill(alleleVariance, 0);
                            //Arrays.fill(averageHetFrequencies, 0);

                            for (int allele = 0; allele < numAlleles; allele++) {
                                Arrays.fill(alleleFrequencies[allele], 0);
                                Arrays.fill(hetFrequencies[allele], 0);
                            }

                            int popIndex = 0;
                            for (Document populationResult : populations) {
                                int sampleSize = populationResult.getInteger(FST_RES_SAMPLESIZE);
                                List<Document> alleles = populationResult.getList(FST_RES_ALLELES, Document.class);

                                for (Document alleleResult : alleles) {
                                    int allele = alleleResult.getInteger(FST_RES_ALLELEID);
                                    alleleFrequencies[allele][popIndex] = alleleResult.getDouble(FST_RES_ALLELEFREQUENCY);
                                    hetFrequencies[allele][popIndex] = alleleResult.getDouble(FST_RES_HETEROZYGOUSFREQUENCY);
                                }

                                sampleSizes[popIndex] = sampleSize;
                                popIndex += 1;
                            }

                            double averageSampleSize = (double)IntStream.of(sampleSizes).sum() / numPopulations;  // n¯ : Average sample size
                            double totalSize = averageSampleSize * numPopulations;  // r × n¯
                            double sampleSizeCorrection = (totalSize - IntStream.of(sampleSizes).mapToDouble(size -> size*size / totalSize).sum() / (numPopulations - 1));  // n_c

                            for (int allele = 0; allele < numAlleles; allele++) {
                                // Compute weighted averages of allele frequencies (p¯) and heterozygote proportions (h¯)
                                double averageAlleleFrequency = 0.0;
                                double averageHetFrequency = 0.0;
                                for (popIndex = 0; popIndex < numPopulations; popIndex++) {
                                    averageAlleleFrequency += sampleSizes[popIndex] * alleleFrequencies[allele][popIndex] / totalSize;
                                    averageHetFrequency += sampleSizes[popIndex] * hetFrequencies[allele][popIndex] / totalSize;
                                }
                                averageAlleleFrequency = new BigDecimal(averageAlleleFrequency).setScale(15, RoundingMode.CEILING).doubleValue();   // round it up because it may look like 0.999999999 but actually mean 1

                                // Compute allele frequency variance (s²)
                                double alleleVariance = 0.0;
                                for (popIndex = 0; popIndex < numPopulations; popIndex++) {
                                    alleleVariance += sampleSizes[popIndex] * Math.pow(alleleFrequencies[allele][popIndex] - averageAlleleFrequency, 2) / (averageSampleSize * (numPopulations - 1));
                                }

                                // a = (n¯/nc) × (s² - (1 / (n¯-1))(p¯(1-p¯) - s²(r-1) / r - h¯/4))
                                double populationVariance = (
                                        (averageSampleSize / sampleSizeCorrection) * (
                                            alleleVariance - (1 / (averageSampleSize - 1)) * (
                                                averageAlleleFrequency * (1 - averageAlleleFrequency) -
                                                (numPopulations - 1) * alleleVariance / numPopulations -
                                                averageHetFrequency / 4
                                            )
                                        )
                                    );

                                // b = (n¯ / (n¯-1))(p¯(1-p¯) - s²(r-1) / r - h¯(2n¯-1)/4n¯)
                                double individualVariance = (averageSampleSize / (averageSampleSize - 1)) * (
                                        averageAlleleFrequency * (1 - averageAlleleFrequency) -
                                        alleleVariance * (numPopulations-1) / numPopulations -
                                        averageHetFrequency * (2*averageSampleSize - 1) / (4 * averageSampleSize)
                                    );

                                // c = h¯/2
                                double gameteVariance = averageHetFrequency / 2;

                                if (!Double.isNaN(populationVariance) && !Double.isNaN(individualVariance) && !Double.isNaN(gameteVariance)) {
                                    weightedFstSum += populationVariance;
                                    fstWeight += populationVariance + individualVariance + gameteVariance;
                                }
                            }
                        }

                        result.put(intervalEntry.getKey(), weightedFstSum / fstWeight);
                        finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                        nTotalTreatedVariantCount.addAndGet(partialCount);
                    }
                    catch (MongoCommandException mce) {
                        progress.setError("Error during Fst calculation: " + mce.getErrorMessage());
                    }
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId()); // important to be sure that all tasks in the group are executed before the queue purges it
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionFst treated " + nTotalTreatedVariantCount.get() + " variants on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");

        progress.markAsComplete();

        return new TreeMap<Long, Double>(result);
    }

    public List<Map<Long, Double>> selectionTajimaD(MgdbChartRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Tajima's D on sequence " + gdr.getDisplayedSequence()});
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        Collection<BasicDBList> variantDataQueries = nTempVarCount == 0 ? varQueryWrapper.getVariantRunDataQueries() : varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();
        final String vrdCollName = mongoTemplate.getCollectionName(VariantRunData.class);
        final boolean useTempColl = (nTempVarCount != 0);
        final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : vrdCollName;
        final ConcurrentHashMap<Long, Double> tajimaD = new ConcurrentHashMap<Long, Double>();
        final ConcurrentHashMap<Long, Double> segregatingSites = new ConcurrentHashMap<Long, Double>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, useTempColl ? usedVarCollName : null)) {
                progress.setError("selectionTajimaD: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filters).");
                return new ArrayList<>();
            }

        List<BasicDBObject> baseQuery = buildTajimaDQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
                segregatingSites.put(intervalEntry.getKey(), 0.0);
                tajimaD.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

        	List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
	            intervalEntry.getValue().append(VariantData.FIELDNAME_KNOWN_ALLELES + ".2", new BasicDBObject("$exists", false)); // exclude multi-allelic variants from MAF calculation
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;
                    Document chunk = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(true).first();  // There's only one interval per query
                    if (chunk != null) {
                        double value = chunk.getDouble(TJD_RES_TAJIMAD);
                        double sites = (double)chunk.getInteger(TJD_RES_SEGREGATINGSITES);
                        segregatingSites.put(intervalEntry.getKey(), sites);
                        tajimaD.put(intervalEntry.getKey(), value);
                    } else {
                        segregatingSites.put(intervalEntry.getKey(), 0.0);
                        tajimaD.put(intervalEntry.getKey(), Double.NaN);
                    }
                    progress.setCurrentStepProgress((short) segregatingSites.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId()); // important to be sure that all tasks in the group are executed before the queue purges it
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionTajimaD treated " + threadsToWaitFor.size() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");

        progress.markAsComplete();

        return Arrays.asList(new TreeMap<>(tajimaD), new TreeMap<>(segregatingSites));
    }

    public Map<Long, Float> selectionMaf(MgdbChartRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
                .orElseGet(() -> {
                    ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "MAF on sequence " + gdr.getDisplayedSequence()});
                    ProgressIndicator.registerProgressIndicator(pi);
                    return pi;
                });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0) {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        Collection<BasicDBList> variantDataQueries = nTempVarCount == 0 ? varQueryWrapper.getVariantRunDataQueries() : varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();
        final String vrdCollName = mongoTemplate.getCollectionName(VariantRunData.class);
        final boolean useTempColl = (nTempVarCount != 0);
        final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : vrdCollName;
        final Map<Long, Float> result = new ConcurrentHashMap<>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, useTempColl ? usedVarCollName : null)) {
                progress.setError("selectionMaf: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filters).");
                return result;
            }

        List<BasicDBObject> baseQuery = buildMafQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
            if (intervalEntry.getValue() == null) {
                result.put(intervalEntry.getKey(), Float.NaN);
                continue;
            }

            List<BasicDBObject> windowQuery = new ArrayList<>(baseQuery);
            if (!useTempColl)
                intervalEntry.getValue().append(VariantData.FIELDNAME_KNOWN_ALLELES + ".2", new BasicDBObject("$exists", false));
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    Document chunk = mongoTemplate.getCollection(usedVarCollName)
                            .aggregate(windowQuery).allowDiskUse(true).first();

                    // Pipeline outputs { _id: null, maf: <value> } — one document per interval
                    if (chunk == null || chunk.get("maf") == null)
                        result.put(intervalEntry.getKey(), Float.NaN);
                    else
                        result.put(intervalEntry.getKey(), 100*chunk.getDouble("maf").floatValue());

                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId());
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionMaf treated " + threadsToWaitFor.size() + " intervals for sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return result;
    }

    private static final String GENOTYPE_DATA_S2_DATA = "dt";
    private static final String GENOTYPE_DATA_S5_SPKEYVAL = "sp";
    private static final String GENOTYPE_DATA_S8_BIO_ENTITY = "be";
    private static final String GENOTYPE_DATA_S7_CALLSET_ID = "cs";
    private static final String GENOTYPE_DATA_S7_GENOTYPE = "gy";
    private static final String GENOTYPE_DATA_S7_POSITION = "ss";
    private static final String GENOTYPE_DATA_S10_VARIANTID = "vi";
    private static final String GENOTYPE_DATA_S10_INDIVIDUALID = "ii";
    private List<BasicDBObject> buildGenotypeDataQuery2(
            MgdbChartRequest gdr,
            boolean useTempColl,
            Map<String, List<Callset>> bioEntityToCallSetListMap,
            boolean keepPosition,
            boolean fGotMultiCallSetIndividuals,
            boolean workWithSamples) throws Exception {

        String[] info = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Integer[] projIDs = Arrays.stream(info[1].split(","))
                .map(Integer::parseInt).toArray(Integer[]::new);

        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        List<BasicDBObject> pipeline = new ArrayList<>();

        // ── STAGE 1 ───────────────────────────────────────────────────────────────
        pipeline.add(null);

        // ── BUILD INDIVIDUAL → GROUP MAPPING ─────────────────────────────────────
        // Derive from gdr.getAllCallSetIds() where index i = group i
        Map<String, Integer> individualToGroupMap = new HashMap<>();
        List<List<String>> allCallSetIds = gdr.getAllCallSetIds();
        for (int groupIdx = 0; groupIdx < allCallSetIds.size(); groupIdx++) {
            for (String callSetId : allCallSetIds.get(groupIdx)) {
                String individualId = callSetId.substring(1 + callSetId.lastIndexOf(Helper.ID_SEPARATOR));
                individualToGroupMap.put(individualId, groupIdx);
            }
        }
        // Default any individual not explicitly mapped to group 0
        for (String ind : bioEntityToCallSetListMap.keySet())
            individualToGroupMap.putIfAbsent(ind, 0);

        // ── BUILD INDIVIDUAL_MAPPING[p][r][i] → {ind, group} ─────────────────────
        Map<Integer, Set<String>> involvedProjectRuns = bioEntityToCallSetListMap.values().stream()
                .flatMap(List::stream)
                .collect(Collectors.groupingBy(Callset::getProjectId,
                        Collectors.mapping(Callset::getRun, Collectors.toSet())));

        // Reverse lookup: "projectId_runName_sampleIndex" → {ind, group}
        Map<String, BasicDBObject> sampleKeyToIndInfo = new HashMap<>();
        for (Map.Entry<String, List<Callset>> entry : bioEntityToCallSetListMap.entrySet()) {
            String individualId = entry.getKey();
            int group = individualToGroupMap.getOrDefault(individualId, 0);
            for (Callset cs : entry.getValue()) {
                String key = cs.getProjectId() + "_" + cs.getRun() + "_" + cs.getCallSetIndexInsideArray();
                sampleKeyToIndInfo.put(key, new BasicDBObject("ind", individualId).append("group", group));
            }
        }

        int spSize = Arrays.stream(projIDs).max(Integer::compareTo).orElse(0) + 1;
        List<Object> individualMapping = new ArrayList<>(Collections.nCopies(spSize, null));

        for (int projId : projIDs) {
            Set<String> involvedRuns = involvedProjectRuns.getOrDefault(projId, Collections.emptySet());

            List<String> orderedRuns = mongoTemplate.findDistinct(
                    new Query(Criteria.where("_id").is(projId)),
                    GenotypingProject.FIELDNAME_RUNS,
                    GenotypingProject.class, String.class);

            List<Object> runsArray = new ArrayList<>(Collections.nCopies(orderedRuns.size(), null));

            for (int rIdx = 0; rIdx < orderedRuns.size(); rIdx++) {
                String runName = orderedRuns.get(rIdx);
                if (!involvedRuns.contains(runName))
                    continue;

                // Size the array for this run from max callsetIndexInsideArray among involved callsets
                int maxIdx = bioEntityToCallSetListMap.values().stream()
                        .flatMap(List::stream)
                        .filter(cs -> cs.getProjectId() == projId && cs.getRun().equals(runName))
                        .mapToInt(Callset::getCallSetIndexInsideArray)
                        .max().orElse(-1);

                if (maxIdx < 0)
                    continue;

                List<Object> samplesArray = new ArrayList<>(Collections.nCopies(maxIdx + 1, null));
                for (int i = 0; i <= maxIdx; i++)
                    samplesArray.set(i, sampleKeyToIndInfo.get(projId + "_" + runName + "_" + i));

                runsArray.set(rIdx, samplesArray);
            }
            individualMapping.set(projId, runsArray);
        }

        // ── STAGE 2 ───────────────────────────────────────────────────────────────
        // Flatten sp[p][r][i] → array of {p, r, i, gt}, skipping nulls at every level
        // ── STAGE 2 ───────────────────────────────────────────────────────────────

        // Innermost: map each index i in gtArr to {p, r, i, gt}
        BasicDBObject gtEntry = new BasicDBObject()
                .append("p", "$$pIdx")
                .append("r", "$$rIdx")
                .append("i", "$$i")
                .append("gt", new BasicDBObject("$arrayElemAt", Arrays.asList("$$gtArr", "$$i")));

        BasicDBObject mapOverGtArr = new BasicDBObject("$map", new BasicDBObject()
                .append("input", new BasicDBObject("$range", Arrays.asList(0, new BasicDBObject("$size", "$$gtArr"))))
                .append("as", "i")
                .append("in", gtEntry));

        // Filter out null gt entries
        BasicDBObject filterNullGt = new BasicDBObject("$filter", new BasicDBObject()
                .append("input", mapOverGtArr)
                .append("as", "e")
                .append("cond", new BasicDBObject("$ne", Arrays.asList("$$e.gt", null))));

        // Inner reduce: iterate over runs in a project
        BasicDBObject gtArrNullCheck = new BasicDBObject("$or", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$$gtArr", null)),
                new BasicDBObject("$not", Arrays.asList(new BasicDBObject("$isArray", "$$gtArr")))));

        BasicDBObject innerReduceIn = new BasicDBObject("$let", new BasicDBObject()
                .append("vars", new BasicDBObject()
                        .append("rIdx", "$$this")
                        .append("gtArr", new BasicDBObject("$arrayElemAt", Arrays.asList("$$runs", "$$this"))))
                .append("in", new BasicDBObject("$cond", new BasicDBObject()
                        .append("if", gtArrNullCheck)
                        .append("then", "$$value")
                        .append("else", new BasicDBObject("$concatArrays", Arrays.asList("$$value", filterNullGt))))));

        BasicDBObject innerReduce = new BasicDBObject("$reduce", new BasicDBObject()
                .append("input", new BasicDBObject("$range", Arrays.asList(0, new BasicDBObject("$size", "$$runs"))))
                .append("initialValue", new ArrayList<>())
                .append("in", innerReduceIn));

        // Outer reduce: iterate over included project indices
        BasicDBObject runsNullCheck = new BasicDBObject("$or", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$$runs", null)),
                new BasicDBObject("$not", Arrays.asList(new BasicDBObject("$isArray", "$$runs")))));

        BasicDBObject outerReduceIn = new BasicDBObject("$let", new BasicDBObject()
                .append("vars", new BasicDBObject()
                        .append("pIdx", "$$this")
                        .append("runs", new BasicDBObject("$arrayElemAt",
                                Arrays.asList("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, "$$this"))))
                .append("in", new BasicDBObject("$cond", new BasicDBObject()
                        .append("if", runsNullCheck)
                        .append("then", "$$value")
                        .append("else", new BasicDBObject("$concatArrays", Arrays.asList("$$value", innerReduce))))));

        BasicDBObject outerReduce = new BasicDBObject("$reduce", new BasicDBObject()
                .append("input", Arrays.asList(projIDs))
                .append("initialValue", new ArrayList<>())
                .append("in", outerReduceIn));

        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject("_gtFlat", outerReduce)));

        // Drop sp immediately — reduces document size before $unwind explosion
        pipeline.add(new BasicDBObject("$unset", VariantRunData.FIELDNAME_SAMPLEGENOTYPES));

        // ── STAGE 3 ───────────────────────────────────────────────────────────────
        pipeline.add(new BasicDBObject("$unwind", "$_gtFlat"));

        // ── STAGE 4 ───────────────────────────────────────────────────────────────

        // Step 1: look up the project array from the mapping
        BasicDBObject lookupProject = new BasicDBObject("$arrayElemAt",
                Arrays.asList(individualMapping, "$_gtFlat.p"));

        // Step 2: look up the run array from the project array
        BasicDBObject lookupRun = new BasicDBObject("$arrayElemAt",
                Arrays.asList(lookupProject, "$_gtFlat.r"));

        // Step 3: look up the {ind, group} object from the run array
        BasicDBObject lookupSample = new BasicDBObject("$arrayElemAt",
                Arrays.asList(lookupRun, "$_gtFlat.i"));

        // Step 4: extract ind and group from the looked-up object
        BasicDBObject indGroupFields = new BasicDBObject()
                .append("ind", "$$info.ind")
                .append("group", "$$info.group");

        // Step 5: wrap in $let to bind the lookup result to $$info
        BasicDBObject letExpr = new BasicDBObject("$let", new BasicDBObject()
                .append("vars", new BasicDBObject("info", lookupSample))
                .append("in", indGroupFields));

        // Step 6: merge into existing _gtFlat document
        BasicDBObject mergeWithIndInfo = new BasicDBObject("$mergeObjects",
                Arrays.asList("$_gtFlat", letExpr));

        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("_gtFlat", mergeWithIndInfo)));

        // ── STAGE 5 ───────────────────────────────────────────────────────────────
        pipeline.add(new BasicDBObject("$match",
                new BasicDBObject("_gtFlat.ind", new BasicDBObject("$ne", null))));

        // ── STAGE 6 ───────────────────────────────────────────────────────────────
        // Count occurrences of each (variant, individual, genotype)
        // p carried with $first — not in _id
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", new BasicDBObject()
                        .append("vi", "$_id." + VariantRunDataId.FIELDNAME_VARIANT_ID)
                        .append("ind", "$_gtFlat.ind")
                        .append("gt", "$_gtFlat.gt")
                        .append("group", "$_gtFlat.group"))
                .append("count", new BasicDBObject("$sum", 1))
                .append("p", new BasicDBObject("$first", "$p"))));

        // ── STAGE 7 ───────────────────────────────────────────────────────────────
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", new BasicDBObject()
                        .append("vi", "$_id.vi")
                        .append("ind", "$_id.ind")
                        .append("group", "$_id.group"))
                .append("gtCounts", new BasicDBObject("$push", new BasicDBObject()
                        .append("gt", "$_id.gt")
                        .append("count", "$count")))
                .append("p", new BasicDBObject("$first", "$p"))));

        // ── STAGE 8 ───────────────────────────────────────────────────────────────
        // Major genotype: unique max → gt, tie → null
        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("_maxCount",
                        new BasicDBObject("$max", "$gtCounts.count"))));

        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("_topGts", new BasicDBObject("$filter", new BasicDBObject()
                        .append("input", "$gtCounts")
                        .append("as", "gc")
                        .append("cond", new BasicDBObject("$eq",
                                Arrays.asList("$$gc.count", "$_maxCount")))))));

        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("majorGt", new BasicDBObject("$cond", new BasicDBObject()
                        .append("if", new BasicDBObject("$eq",
                                Arrays.asList(new BasicDBObject("$size", "$_topGts"), 1)))
                        .append("then", new BasicDBObject("$arrayElemAt",
                                Arrays.asList("$_topGts.gt", 0)))
                        .append("else", null)))));

        // Drop tied individuals — nothing in $genotypes will ever be null downstream
        pipeline.add(new BasicDBObject("$match",
                new BasicDBObject("majorGt", new BasicDBObject("$ne", null))));

        // ── STAGE 9 ───────────────────────────────────────────────────────────────
        // Collect major genotypes per (variant, group)
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", new BasicDBObject()
                        .append("vi", "$_id.vi")
                        .append("group", "$_id.group"))
                .append("genotypes", new BasicDBObject("$push", "$majorGt"))
                .append("p", new BasicDBObject("$first", "$p"))));

        // ── STAGE 10 ──────────────────────────────────────────────────────────────
        // Assemble all groups into one document per variant
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", "$_id.vi")
                .append("groups", new BasicDBObject("$push", new BasicDBObject()
                        .append("group", "$_id.group")
                        .append("genotypes", "$genotypes")))
                .append("p", new BasicDBObject("$first", "$p"))));

        // ── STAGE 11 ──────────────────────────────────────────────────────────────
        BasicDBObject projectFields = new BasicDBObject("_id", 1).append("groups", 1);
        if (keepPosition)
            projectFields.append("p", 1);
        pipeline.add(new BasicDBObject("$project", projectFields));

        LOG.debug("new model genotype query: " + pipeline);
        return pipeline;
    }
    private List<BasicDBObject> buildGenotypeDataQuery(MgdbChartRequest gdr, boolean useTempColl, Map<String, List<Callset>> bioEntityToCallSetListMap, boolean keepPosition, boolean fGotMultiCallSetIndividuals, boolean workWithSamples) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Integer[] projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toArray(Integer[]::new);
        boolean fWillNeedToMergeObjects = projIDs.length > 1;   // if multiple projects or runs, we will need to merge objects to have just one record per variant
        if (!fWillNeedToMergeObjects && bioEntityToCallSetListMap.values().stream().flatMap(List::stream).map(cs -> cs.getRun()).distinct().count() > 1)
            fWillNeedToMergeObjects = true; // even in a single project, multiple runs involved

        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        List<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();

        // Stage 1 : placeholder for initial match stage
        pipeline.add(null);

        if (useTempColl) {
            // Stage 2 : Lookup from temp collection to variantRunData
            BasicDBObject lookup = new BasicDBObject();
            lookup.put("from", mongoTemplate.getCollectionName(VariantRunData.class));
            lookup.put("localField", "_id");
            lookup.put("foreignField", "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
            lookup.put("as", GENOTYPE_DATA_S2_DATA);
            pipeline.add(new BasicDBObject("$lookup", lookup));

            // Stage 3 : Unwind data
            pipeline.add(new BasicDBObject("$unwind", "$" + GENOTYPE_DATA_S2_DATA));
        }

        // Stage 4 : Keep only the right projects / runs
        ArrayList<BasicDBObject> projectFilterList = new ArrayList<>();
        Map<Integer, Set<String>> involvedProjectRuns = bioEntityToCallSetListMap.values().stream().flatMap(List::stream).collect(Collectors.groupingBy(Callset::getProjectId,Collectors.mapping(Callset::getRun, Collectors.toSet())));
        for (int projId : involvedProjectRuns.keySet()) {
        	Set<String> projectInvolvedRuns = involvedProjectRuns.get(projId);
            BasicDBObject projectFilter = new BasicDBObject((useTempColl ? GENOTYPE_DATA_S2_DATA + "." : "") + "_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, projId);
            // if not all of this project's runs are involved, only match the required ones
            boolean fNotAllRunsNeeded = projectInvolvedRuns.size() != mongoTemplate.findDistinct(new Query(Criteria.where("_id").is(projId)), GenotypingProject.FIELDNAME_RUNS, GenotypingProject.class, String.class).size();
            projectFilterList.add(fNotAllRunsNeeded ? new BasicDBObject("$and", Arrays.asList(projectFilter, new BasicDBObject((useTempColl ? GENOTYPE_DATA_S2_DATA + "." : "") + "_id." + VariantRunDataId.FIELDNAME_RUNNAME, new BasicDBObject("$in", projectInvolvedRuns)))) : projectFilter);
        }

        if (!projectFilterList.isEmpty())
        	pipeline.add(new BasicDBObject("$match", projectFilterList.size() == 1 ? projectFilterList.get(0) : new BasicDBObject("$or", projectFilterList)));

        if (fGotMultiCallSetIndividuals) {
            // Stage 5 : Convert samples to an array
            pipeline.add(new BasicDBObject("$set", new BasicDBObject(GENOTYPE_DATA_S5_SPKEYVAL, new BasicDBObject("$objectToArray", "$" + (useTempColl ? GENOTYPE_DATA_S2_DATA + "." + VariantRunData.FIELDNAME_SAMPLEGENOTYPES : VariantRunData.FIELDNAME_SAMPLEGENOTYPES)))));

            // Stage 6 : Unwind samples
            pipeline.add(new BasicDBObject("$unwind", "$" + GENOTYPE_DATA_S5_SPKEYVAL));

            // Stage 7 : Convert key and value
            BasicDBObject spkeyval = new BasicDBObject();
            spkeyval.put(GENOTYPE_DATA_S7_CALLSET_ID, new BasicDBObject("$toInt", "$" + GENOTYPE_DATA_S5_SPKEYVAL + ".k"));
            spkeyval.put(GENOTYPE_DATA_S7_GENOTYPE, "$" + GENOTYPE_DATA_S5_SPKEYVAL + ".v.gt");
            if (keepPosition)
                spkeyval.put(GENOTYPE_DATA_S7_POSITION, "$" + Assembly.getThreadBoundVariantRefPosPath() + ".ss");
            pipeline.add(new BasicDBObject("$project", spkeyval));

            // Stage 8 : Keep only the callsets we are interested in
            List<Integer> allCallsetIDs = mongoTemplate.findDistinct(new Query(), GenotypingSample.FIELDNAME_CALLSETS + "._id", GenotypingSample.class, Integer.class);
            List<Integer> involvedCallsetIDs = bioEntityToCallSetListMap.values().stream().flatMap(Collection::stream).map(cs -> cs.getId()).distinct().toList();
            if (allCallsetIDs.size() != involvedCallsetIDs.size())	// if all are involved don't bother filtering
            	pipeline.add(new BasicDBObject("$match", new BasicDBObject(GENOTYPE_DATA_S7_CALLSET_ID, involvedCallsetIDs.size() < allCallsetIDs.size() / 2 ? new BasicDBObject("$in", involvedCallsetIDs) : new BasicDBObject("$not", new BasicDBObject("$in", CollectionUtils.disjunction(allCallsetIDs, involvedCallsetIDs))))));

            // Stage 9 : Create bio-entity field using $switch
            BasicDBList branches = new BasicDBList();
            for (Map.Entry<String, List<Callset>> entry : bioEntityToCallSetListMap.entrySet()) {
                if (entry.getValue().size() > 1)
                    branches.add(new BasicDBObject().append("case", new BasicDBObject("$in", new Object[]{"$" + GENOTYPE_DATA_S7_CALLSET_ID, entry.getValue().stream().map(cs -> cs.getId()).toList()})).append("then", entry.getKey()));
            }
            BasicDBObject switchExpr = new BasicDBObject().append("branches", branches).append("default", "$" + GENOTYPE_DATA_S7_CALLSET_ID);
            pipeline.add(new BasicDBObject("$addFields",new BasicDBObject(GENOTYPE_DATA_S8_BIO_ENTITY , new BasicDBObject("$switch", switchExpr))));

            // Stage 10 : Regroup individual runs
            BasicDBObject individualGroup = new BasicDBObject();
            BasicDBObject individualGroupId = new BasicDBObject();
            individualGroupId.put(GENOTYPE_DATA_S10_VARIANTID, "$_id" + (!useTempColl? "." + FST_S22_VARIANTID : ""));
            individualGroupId.put(GENOTYPE_DATA_S10_INDIVIDUALID, "$" + GENOTYPE_DATA_S8_BIO_ENTITY );
            individualGroup.put("_id", individualGroupId);
            individualGroup.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$addToSet", "$" + GENOTYPE_DATA_S7_GENOTYPE));
            individualGroup.put(VariantRunDataId.FIELDNAME_PROJECT_ID, new BasicDBObject("$addToSet", "$_id." + VariantRunDataId.FIELDNAME_PROJECT_ID));
            if (keepPosition)
                individualGroup.put(GENOTYPE_DATA_S7_POSITION, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S7_POSITION));
            pipeline.add(new BasicDBObject("$group", individualGroup));

            // Stage 11 : Weed out incoherent genotypes
            pipeline.add(new BasicDBObject("$match", new BasicDBObject(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$size", 1))));

            // Stage 12 : Convert back to sp object
            pipeline.add(new BasicDBObject("$project", new BasicDBObject("_id", "$_id." + FST_S22_VARIANTID)
            		.append(VariantRunDataId.FIELDNAME_PROJECT_ID, 1).append(GENOTYPE_DATA_S7_POSITION, 1)
            		.append(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$arrayToObject", Arrays.asList(Arrays.asList(new BasicDBObject()
																        .append("k", new BasicDBObject("$toString", "$_id." + GENOTYPE_DATA_S10_INDIVIDUALID))
																        .append("v", new BasicDBObject(TJD_S18_GENOTYPE, new BasicDBObject("$arrayElemAt", Arrays.asList("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, 0))))))))
            ));
        }

        if (fWillNeedToMergeObjects || useTempColl) {  // merge genotypes into a single document per variant
            BasicDBObject projectMergeGroup = new BasicDBObject();
            projectMergeGroup.put("_id", "$_id" + (!useTempColl && !fGotMultiCallSetIndividuals ? "." + FST_S22_VARIANTID : ""));
            projectMergeGroup.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$mergeObjects", "$" + (useTempColl && !fGotMultiCallSetIndividuals ? GENOTYPE_DATA_S2_DATA + "." : "") + VariantRunData.FIELDNAME_SAMPLEGENOTYPES));
            pipeline.add(new BasicDBObject("$group", projectMergeGroup));
        }

        LOG.debug("base query " + pipeline);
        return pipeline;
    }

    private static final String FST_S14_POPULATIONGENOTYPES = "pg";
    private static final String FST_S15_POPULATION = "pp";
    private static final String FST_S19_GENOTYPE = "gn";
    private static final String FST_S20_HETEROZYGOTE = "ht";
    private static final String FST_S22_VARIANTID = "vi";
    private static final String FST_S22_POPULATIONID = "pi";
    private static final String FST_S22_ALLELECOUNT = "ac";
    private static final String FST_S22_HETEROZYGOUSCOUNT = "hc";

    private static final String FST_RES_SAMPLESIZE = "ss";
    private static final String FST_RES_ALLELEID = "al";
    private static final String FST_RES_ALLELEMAX = "am";
    private static final String FST_RES_ALLELEFREQUENCY = "af";
    private static final String FST_RES_HETEROZYGOUSFREQUENCY = "hf";
    private static final String FST_RES_ALLELES = "as";
    private static final String FST_RES_POPULATIONS = "ps";

    private List<BasicDBObject> buildFstQuery(MgdbChartRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        LOG.debug("Calculating Fst on " + info[0] + " with group sizes: " + gdr.getAllCallSetIds().stream().map(t -> t.size()).toList());

        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        // ── unchanged from original ───────────────────────────────────────────────
        List<Collection<String>> selectedMaterial = new ArrayList<>();
        List<List<String>> ga4ghCallsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < ga4ghCallsetIds.size(); i++)
            selectedMaterial.add(ga4ghCallsetIds.get(i).isEmpty()
                    ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs))
                    : ga4ghCallsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        for (Collection<String> group : selectedMaterial)
            if (!workWithSamples)
                individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, group));
            else
                individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, group));

        boolean fGotMultiCallSetIndividuals = individualOrSampleToCallSetListMap.values().stream()
                .anyMatch(spList -> spList.size() > 1);
        // ─────────────────────────────────────────────────────────────────────────

        List<BasicDBObject> pipeline = buildGenotypeDataQuery2(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // ── FST STAGE 1 ───────────────────────────────────────────────────────────
        // One document per group per variant
        pipeline.add(new BasicDBObject("$unwind", "$groups"));

        // ── FST STAGE 2 ───────────────────────────────────────────────────────────
        // Compute raw counts directly from integer encoding:
        // abs(gt) mod 100 = alt allele count
        // floor(abs(gt) / 100) = ploidy (total alleles per individual)
        // gt < 0 = heterozygote
        BasicDBObject mapAltCount = new BasicDBObject("$map", new BasicDBObject()
                .append("input", "$groups.genotypes")
                .append("as", "g")
                .append("in", new BasicDBObject("$mod", Arrays.asList(
                        new BasicDBObject("$abs", "$$g"), 100))));

        BasicDBObject mapPloidy = new BasicDBObject("$map", new BasicDBObject()
                .append("input", "$groups.genotypes")
                .append("as", "g")
                .append("in", new BasicDBObject("$floor", new BasicDBObject("$divide", Arrays.asList(
                        new BasicDBObject("$abs", "$$g"), 100)))));

        BasicDBObject mapIsHet = new BasicDBObject("$map", new BasicDBObject()
                .append("input", "$groups.genotypes")
                .append("as", "g")
                .append("in", new BasicDBObject("$cond", Arrays.asList(
                        new BasicDBObject("$lt", Arrays.asList("$$g", 0)), 1, 0))));

        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject("_fstRaw", new BasicDBObject()
                .append(FST_RES_SAMPLESIZE, new BasicDBObject("$size", "$groups.genotypes"))
                .append("ta", new BasicDBObject("$sum", mapPloidy))     // total alleles
                .append("aa", new BasicDBObject("$sum", mapAltCount))   // alt allele count
                .append("hc", new BasicDBObject("$sum", mapIsHet)))));  // heterozygote count

        // ── FST STAGE 3 ───────────────────────────────────────────────────────────
        // Compute allele and het frequencies, build population document
        // Matches the structure the Java Fst algorithm already reads:
        // { ss: sampleSize, as: [{al:0, af:refFreq, hf:hetFreq}, {al:1, af:altFreq, hf:hetFreq}] }
        BasicDBObject refFreq = new BasicDBObject("$cond", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$_fstRaw.ta", 0)), 0,
                new BasicDBObject("$divide", Arrays.asList(
                        new BasicDBObject("$subtract", Arrays.asList("$_fstRaw.ta", "$_fstRaw.aa")),
                        "$_fstRaw.ta"))));

        BasicDBObject altFreq = new BasicDBObject("$cond", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$_fstRaw.ta", 0)), 0,
                new BasicDBObject("$divide", Arrays.asList("$_fstRaw.aa", "$_fstRaw.ta"))));

        BasicDBObject hetFreq = new BasicDBObject("$cond", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$_fstRaw." + FST_RES_SAMPLESIZE, 0)), 0,
                new BasicDBObject("$divide", Arrays.asList(
                        "$_fstRaw.hc", "$_fstRaw." + FST_RES_SAMPLESIZE))));

        BasicDBObject allele0 = new BasicDBObject()
                .append(FST_RES_ALLELEID, 0)
                .append(FST_RES_ALLELEFREQUENCY, refFreq)
                .append(FST_RES_HETEROZYGOUSFREQUENCY, hetFreq);

        BasicDBObject allele1 = new BasicDBObject()
                .append(FST_RES_ALLELEID, 1)
                .append(FST_RES_ALLELEFREQUENCY, altFreq)
                .append(FST_RES_HETEROZYGOUSFREQUENCY, hetFreq);

        BasicDBObject population = new BasicDBObject()
                .append(FST_RES_SAMPLESIZE, "$_fstRaw." + FST_RES_SAMPLESIZE)
                .append(FST_RES_ALLELES, Arrays.asList(allele0, allele1));

        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject("_pop", population)));

        // ── FST STAGE 4 ───────────────────────────────────────────────────────────
        // Group by variant, collect all population documents.
        // alleleMax = 1 (biallelic already enforced by Stage 1 match filter)
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", "$_id")
                .append(FST_RES_POPULATIONS, new BasicDBObject("$push", "$_pop"))
                .append(FST_RES_ALLELEMAX, new BasicDBObject("$first", 1))));

        LOG.debug("Fst query " + pipeline);
        return pipeline;
    }

    private static final String TJD_S14_GENOTYPES = "gl";
    private static final String TJD_S15_SAMPLESIZE = "sz";
    private static final String TJD_S18_GENOTYPE = SampleGenotype.FIELDNAME_GENOTYPECODE;
    private static final String TJD_S20_VARIANTID = "vi";
    private static final String TJD_S20_ALLELEID = "ai";
    private static final String TJD_S20_ALLELECOUNT = "ac";
    private static final String TJD_S21_NUMALLELES = "na";
    private static final String TJD_S23_ALLELEFREQUENCY = "k";
    private static final String TJD_S24_FREQUENCYSUM = "kc";

    private static final String TJD_RES_SEGREGATINGSITES = "sg";
    private static final String TJD_RES_TAJIMAD = "tjd";

    private List<BasicDBObject> buildTajimaDQuery(MgdbChartRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        // ── unchanged from original ───────────────────────────────────────────────
        HashSet<String> selectedMaterial = new HashSet<String>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));

        final int sampleSize = 2 * selectedMaterial.size();

        double a1 = 0, a2 = 0;
        for (int i = 1; i < sampleSize; i++) {
            a1 += 1.0 / i;
            a2 += 1.0 / (i * i);
        }
        double b1 = (double)(sampleSize + 1) / (double)(3 * (sampleSize - 1));
        double b2 = 2.0 * (sampleSize * sampleSize + sampleSize + 3) / (9.0 * sampleSize * (sampleSize - 1));
        double c1 = b1 - 1 / a1;
        double c2 = b2 - (double)(sampleSize + 2) / (a1 * sampleSize) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        // ─────────────────────────────────────────────────────────────────────────

        boolean fGotMultiCallSetIndividuals = (individualOrSampleToCallSetListMap.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent());
        List<BasicDBObject> pipeline = buildGenotypeDataQuery2(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // ── TAJIMAD STAGE 1 ───────────────────────────────────────────────────────
        // Extract genotypes from first group (single group for Tajima's D)
        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("genotypes",
                        new BasicDBObject("$arrayElemAt", Arrays.asList("$groups.genotypes", 0)))));
        pipeline.add(new BasicDBObject("$unset", "groups"));

        // ── TAJIMAD STAGE 2 ───────────────────────────────────────────────────────
        // Compute alt allele count and per-variant sample size from integer encoding
        // abs(gt) mod 100 = alt allele count per individual
        BasicDBObject mapAltCount = new BasicDBObject("$map", new BasicDBObject()
                .append("input", "$genotypes")
                .append("as", "g")
                .append("in", new BasicDBObject("$mod", Arrays.asList(
                        new BasicDBObject("$abs", "$$g"), 100))));

        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject()
                .append(TJD_S20_ALLELECOUNT, new BasicDBObject("$sum", mapAltCount))
                .append(TJD_S15_SAMPLESIZE, new BasicDBObject("$size", "$genotypes"))));

        // ── TAJIMAD STAGE 3 ───────────────────────────────────────────────────────
        // Filter non-segregating variants (monomorphic ref or monomorphic alt)
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$expr",
                new BasicDBObject("$and", Arrays.asList(
                        new BasicDBObject("$gt", Arrays.asList("$" + TJD_S20_ALLELECOUNT, 0)),
                        new BasicDBObject("$lt", Arrays.asList("$" + TJD_S20_ALLELECOUNT,
                                new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S15_SAMPLESIZE, 2)))))))));

        // ── TAJIMAD STAGE 4 ───────────────────────────────────────────────────────
        // Compute p(1-p) per variant — unchanged formula from original Stage 23
        BasicDBObject alleleFrequency = new BasicDBObject();
        alleleFrequency.put("vars", new BasicDBObject("freq",
                new BasicDBObject("$divide", Arrays.asList(
                        "$" + TJD_S20_ALLELECOUNT,
                        new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S15_SAMPLESIZE, 2))))));
        alleleFrequency.put("in", new BasicDBObject("$multiply", Arrays.asList(
                "$$freq", new BasicDBObject("$subtract", Arrays.asList(1, "$$freq")))));

        BasicDBObject frequencyProjection = new BasicDBObject();
        frequencyProjection.put(TJD_S23_ALLELEFREQUENCY, new BasicDBObject("$let", alleleFrequency));
        pipeline.add(new BasicDBObject("$project", frequencyProjection));

        // ── TAJIMAD STAGE 5 ───────────────────────────────────────────────────────
        // Group interval — unchanged from original Stage 24
        BasicDBObject outputGroup = new BasicDBObject();
        outputGroup.put("_id", 1);
        outputGroup.put(TJD_S24_FREQUENCYSUM, new BasicDBObject("$sum", "$" + TJD_S23_ALLELEFREQUENCY));
        outputGroup.put(TJD_RES_SEGREGATINGSITES, new BasicDBObject("$sum", 1));
        pipeline.add(new BasicDBObject("$group", outputGroup));

        // ── TAJIMAD STAGE 6 ───────────────────────────────────────────────────────
        // Compute Tajima's D — unchanged from original Stage 25
        BasicDBObject finalProject = new BasicDBObject();
        finalProject.put(TJD_RES_SEGREGATINGSITES, 1);
        finalProject.put(TJD_RES_TAJIMAD, new BasicDBObject("$divide", Arrays.asList(
                new BasicDBObject("$subtract", Arrays.asList(
                        new BasicDBObject("$divide", Arrays.asList(
                                new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S24_FREQUENCYSUM, 2 * sampleSize)),
                                sampleSize - 1)),
                        new BasicDBObject("$divide", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, a1)))),
                new BasicDBObject("$sqrt", new BasicDBObject("$abs", new BasicDBObject("$add", Arrays.asList(
                        new BasicDBObject("$multiply", Arrays.asList(e1, "$" + TJD_RES_SEGREGATINGSITES)),
                        new BasicDBObject("$multiply", Arrays.asList(
                                new BasicDBObject("$multiply", Arrays.asList(e2, "$" + TJD_RES_SEGREGATINGSITES)),
                                new BasicDBObject("$subtract", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, 1)))))))))));
        pipeline.add(new BasicDBObject("$project", finalProject));
        LOG.debug("Tajima's D query: " + pipeline);

        return pipeline;
    }
    private List<BasicDBObject> buildMafQuery(MgdbChartRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        String[] info = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(Integer::parseInt).toList();
        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);

        // Build selected material from all groups
        HashSet<String> selectedMaterial = new HashSet<>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty()
                    ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs))
                    : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));

        boolean fGotMultiCallSetIndividuals = individualOrSampleToCallSetListMap.values().stream()
                .anyMatch(spList -> spList.size() > 1);

        // Base pipeline from new model
        List<BasicDBObject> pipeline = buildGenotypeDataQuery2(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // ── MAF STAGE 1 ───────────────────────────────────────────────────────────
        // Extract genotypes from first group (MAF works on a single group)
        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("genotypes",
                        new BasicDBObject("$arrayElemAt", Arrays.asList("$groups.genotypes", 0)))));
        pipeline.add(new BasicDBObject("$unset", "groups"));

        // ── MAF STAGE 2 ───────────────────────────────────────────────────────────
        // Filter out null genotypes (ties or missing) before computing stats
                BasicDBObject filterNulls = new BasicDBObject("$filter", new BasicDBObject()
                .append("input", "$genotypes")
                .append("as", "g")
                .append("cond", new BasicDBObject("$ne", Arrays.asList("$$g", null))));

        // Compute alt allele count and total alleles from genotype encoding
    // encoding: abs(gt) mod 100 = alt allele count, floor(abs(gt) / 100) = ploidy
        BasicDBObject altCalc = new BasicDBObject("$abs",
                new BasicDBObject("$mod", Arrays.asList("$$this", 100)));

        BasicDBObject totalCalc = new BasicDBObject("$floor",
                new BasicDBObject("$divide", Arrays.asList(
                        new BasicDBObject("$abs", "$$this"), 100)));

        BasicDBObject reduceIn = new BasicDBObject()
                .append("alt", new BasicDBObject("$add", Arrays.asList("$$value.alt", altCalc)))
                .append("total", new BasicDBObject("$add", Arrays.asList("$$value.total", totalCalc)));

        BasicDBObject statsReduce = new BasicDBObject("$reduce", new BasicDBObject()
                .append("input", filterNulls)        // ← filtered input, not raw "$genotypes"
                .append("initialValue", new BasicDBObject("alt", 0).append("total", 0))
                .append("in", reduceIn));

        pipeline.add(new BasicDBObject("$project", new BasicDBObject("stats", statsReduce)));

        // ── MAF STAGE 3 ───────────────────────────────────────────────────────────
        // Compute MAF: if freq > 0.5 use 1-freq, else freq. Null if no data.
        BasicDBObject freqVar = new BasicDBObject("$divide",
                Arrays.asList("$stats.alt", "$stats.total"));

        BasicDBObject mafFromFreq = new BasicDBObject("$cond", Arrays.asList(
                new BasicDBObject("$gt", Arrays.asList("$$freq", 0.5)),
                new BasicDBObject("$subtract", Arrays.asList(1, "$$freq")),
                "$$freq"));

        BasicDBObject mafCalc = new BasicDBObject("$cond", Arrays.asList(
                new BasicDBObject("$eq", Arrays.asList("$stats.total", 0)),
                null,
                new BasicDBObject("$let", new BasicDBObject()
                        .append("vars", new BasicDBObject("freq", freqVar))
                        .append("in", mafFromFreq))));

        pipeline.add(new BasicDBObject("$project", new BasicDBObject("maf", mafCalc)));

        // ── MAF STAGE 4 ───────────────────────────────────────────────────────────
        // Average MAF across all variants in the interval
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", null)
                .append("maf", new BasicDBObject("$avg", "$maf"))));

        LOG.debug("MAF query " + pipeline);
        return pipeline;
    }

    public Map<Long, Float> selectionMissingData(MgdbChartRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "missing data percentage on sequence " + gdr.getDisplayedSequence()});
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        Collection<BasicDBList> variantDataQueries = nTempVarCount == 0 ? varQueryWrapper.getVariantRunDataQueries() : varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();
        final String vrdCollName = mongoTemplate.getCollectionName(VariantRunData.class);
        final boolean useTempColl = (nTempVarCount != 0);
        final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : vrdCollName;
        final Map<Long, Float> result = new ConcurrentHashMap<>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, useTempColl ? usedVarCollName : null)) {
                progress.setError("selectionMissingData: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filters).");
                return result;
            }
        
        List<BasicDBObject> baseQuery = buildMissingDataQuery(gdr, useTempColl, workWithSamples);
        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();
        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
            //LOG.debug("missing data window query: " + windowQuery);
            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    long totalVariants = mongoTemplate.getCollection(usedVarCollName)
                            .countDocuments(intervalEntry.getValue());

                    Document chunk = mongoTemplate.getCollection(usedVarCollName)
                            .aggregate(windowQuery).allowDiskUse(true).first();

                    if (totalVariants == 0)
                        result.put(intervalEntry.getKey(), Float.NaN);       // truly empty interval
                    else if (chunk == null)
                        result.put(intervalEntry.getKey(), 100f);             // all variants were 100% missing
                    else {
                        long sumNmc = chunk.getInteger("sumNmc");
                        long totalIndividuals = chunk.getInteger("totalIndividuals");
                        long totalPossible = totalVariants * totalIndividuals;
                        result.put(intervalEntry.getKey(), totalPossible == 0 ? 0f
                                : (float)(totalPossible - sumNmc) / totalPossible * 100);
                    }
                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId());
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionMissingData treated " + threadsToWaitFor.size() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<>(result);
    }

    /**
     * Builds aggregation pipeline to calculate missing data percentage
     * Missing data is defined as when the genotype field is null or empty
     *
     * @param gdr the density request containing query parameters
     * @param useTempColl whether to use a temporary collection
     * @param workWithSamples whether to work with samples or individuals
     * @return pipeline stages for missing data calculation
     * @throws Exception if an error occurs
     */
    private List<BasicDBObject> buildMissingDataQuery(MgdbChartRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        // -------------------- Extract module and project IDs --------------------
        String[] info = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(Integer::parseInt).toList();

        // -------------------- Determine all individuals/samples to consider --------------------
        HashSet<String> selectedMaterial = new HashSet<>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++) {
            List<String> currentCallsetIds = callsetIds.get(i);
            if (currentCallsetIds.isEmpty()) {
                if (workWithSamples) {
                    selectedMaterial.addAll(MgdbDao.getProjectSamples(info[0], projIDs));
                } else {
                    selectedMaterial.addAll(MgdbDao.getProjectIndividuals(info[0], projIDs));
                }
            } else {
                selectedMaterial.addAll(currentCallsetIds.stream()
                        .map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR)))
                        .collect(Collectors.toSet()));
            }
        }
        final int totalIndividuals = selectedMaterial.size(); // known universe

        // -------------------- Build mapping individual/sample -> callset IDs --------------------
        HashMap<String, List<Callset>> individualOrSampleToCallSetListMap = new HashMap<>();
        if (!workWithSamples) {
            individualOrSampleToCallSetListMap.putAll(
                    MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        } else {
            individualOrSampleToCallSetListMap.putAll(
                    MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));
        }

        // -------------------- Build allowed keys for filtering $sp --------------------
        Set<String> allowedKeySet = new HashSet<>();
        for (Map.Entry<String, List<Callset>> entry : individualOrSampleToCallSetListMap.entrySet())
            if (entry.getValue().size() > 1)
            	allowedKeySet.add(entry.getKey());
            else
            	allowedKeySet.add(String.valueOf(entry.getValue().get(0).getId()));

        // -------------------- Start building aggregation pipeline --------------------
        List<BasicDBObject> pipeline = buildGenotypeDataQuery2(gdr, useTempColl, individualOrSampleToCallSetListMap, true,
                individualOrSampleToCallSetListMap.values().stream().anyMatch(spList -> spList.size() > 1),
                workWithSamples);

        // ── MISSING STAGE 1 ───────────────────────────────────────────────────────
        // Extract genotypes from first group
        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("genotypes",
                        new BasicDBObject("$arrayElemAt", Arrays.asList("$groups.genotypes", 0)))));
        pipeline.add(new BasicDBObject("$unset", "groups"));

        // ── MISSING STAGE 2 ───────────────────────────────────────────────────────
        // nmc = number of individuals with data for this variant
        pipeline.add(new BasicDBObject("$project",
                new BasicDBObject("nmc", new BasicDBObject("$size", "$genotypes"))));

        // ── MISSING STAGE 3 ───────────────────────────────────────────────────────
        // Average missing rate: (totalIndividuals - nmc) / totalIndividuals
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", null)
                .append("sumNmc", new BasicDBObject("$sum", "$nmc"))
                .append("variantCount", new BasicDBObject("$sum", 1))
                .append("totalIndividuals", new BasicDBObject("$first", totalIndividuals))));  // ← carry it through

        LOG.debug("Missing data query: " + pipeline);
        return pipeline;
    }

    private static final String HET_S14_POPULATIONGENOTYPES = "pg";
    private static final String HET_S15_POPULATION = "pp";
    private static final String HET_S16_SAMPLESIZE = "ss";
    private static final String HET_S17_GENOTYPE = "gt";
    private static final String HET_S18_ISHETEROZYGOTE = "ht";

    private static final String HET_RES_HETEROZYGOSITY = "he";
    private static final String HET_RES_POPULATIONS = "ps";

    /**
     * Builds aggregation pipeline to calculate heterozygosity rate across populations
     * Heterozygosity is the proportion of heterozygous genotypes per variant
     *
     * @param gdr the density request containing query parameters
     * @param useTempColl whether to use a temporary collection
     * @param workWithSamples whether to work with samples or individuals
     * @return pipeline stages for heterozygosity calculation
     * @throws Exception if an error occurs
     */
    private List<BasicDBObject> buildHeterozygosityQuery(MgdbChartRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        String[] info = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(Integer::parseInt).toList();

        HashSet<String> selectedMaterial = new HashSet<>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty()
                    ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs))
                    : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));

        boolean fGotMultiCallSetIndividuals = individualOrSampleToCallSetListMap.values().stream()
                .anyMatch(spList -> spList.size() > 1);

        List<BasicDBObject> pipeline = buildGenotypeDataQuery2(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // ── HET STAGE 1 ───────────────────────────────────────────────────────────
        // Extract genotypes from first group
        pipeline.add(new BasicDBObject("$addFields",
                new BasicDBObject("genotypes",
                        new BasicDBObject("$arrayElemAt", Arrays.asList("$groups.genotypes", 0)))));
        pipeline.add(new BasicDBObject("$unset", "groups"));

        // ── HET STAGE 2 ───────────────────────────────────────────────────────────
        // Count heterozygotes: negative genotypes in the encoding are heterozygotes
        BasicDBObject filterHeterozygotes = new BasicDBObject("$filter", new BasicDBObject()
                .append("input", "$genotypes")
                .append("cond", new BasicDBObject("$lt", Arrays.asList("$$this", 0))));

        pipeline.add(new BasicDBObject("$project", new BasicDBObject()
                .append("hc", new BasicDBObject("$size", filterHeterozygotes))
                .append("t", new BasicDBObject("$size", "$genotypes"))));

        // ── HET STAGE 3 ───────────────────────────────────────────────────────────
        // Average heterozygote count across all variants in the interval
        pipeline.add(new BasicDBObject("$group", new BasicDBObject()
                .append("_id", null)
                .append("h", new BasicDBObject("$avg",
                        new BasicDBObject("$divide", Arrays.asList("$hc", "$t"))))));

        LOG.debug("Heterozygosity query " + pipeline);
        return pipeline;
    }

    public Map<Long, Float> selectionHeterozygosity(MgdbChartRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
                .orElseGet(() -> {
                    ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "heterozygosity percentage on sequence " + gdr.getDisplayedSequence()});
                    ProgressIndicator.registerProgressIndicator(pi);
                    return pi;
                });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0) {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        Collection<BasicDBList> variantDataQueries = nTempVarCount == 0 ? varQueryWrapper.getVariantRunDataQueries() : varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();
        final boolean useTempColl = (nTempVarCount != 0);
        final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : mongoTemplate.getCollectionName(VariantRunData.class);
        final Map<Long, Float> result = new ConcurrentHashMap<>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, useTempColl ? usedVarCollName : null)) {
                progress.setError("selectionHeterozygosity: Unable to find default position range, either there are no results in the current selection for this sequence, or results are not in sync with interface filters (in which case, try re-applying filters).");
                return result;
            }

        List<BasicDBObject> baseQuery = buildHeterozygosityQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
            if (intervalEntry.getValue() == null) {
                result.put(intervalEntry.getKey(), Float.NaN);
                continue;
            }

            List<BasicDBObject> windowQuery = new ArrayList<>(baseQuery);
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
            LOG.debug("window query: " + windowQuery);
            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    Document chunk = mongoTemplate.getCollection(usedVarCollName)
                            .aggregate(windowQuery).allowDiskUse(true).first();

                    // Pipeline outputs { _id: null, het_avg: <value> } — one document per interval
                    if (chunk == null || chunk.get("h") == null)
                        result.put(intervalEntry.getKey(), Float.NaN);
                    else
                        result.put(intervalEntry.getKey(), 100*chunk.getDouble("h").floatValue());

                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(progress.getProcessId(), t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(progress.getProcessId());
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionHeterozygosity treated " + threadsToWaitFor.size() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<>(result);
    }

    private BasicDBList getFullPathToGenotypes(Collection<String> selectedMaterial, Map<String, List<Callset>> individualOrSampleToCallSetListMap, boolean fGotMultiCallSetIndividuals) {
        BasicDBList result = new BasicDBList();
        Iterator<String> indIt = selectedMaterial.iterator();
        while (indIt.hasNext()) {
            String individualOrSample = indIt.next(), materialRefToUse = null;
            List<Callset> callSets = individualOrSampleToCallSetListMap.get(individualOrSample);
           	materialRefToUse = callSets.size() == 1 ? "" + callSets.get(0).getId() : individualOrSample;
            String pathToGT = materialRefToUse + "." + TJD_S18_GENOTYPE;
            String fullPathToGT = "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES/* + (int) ((individualSample.getId() - 1) / 100)*/ + "." + pathToGT;
            result.add(fullPathToGT);
        }
        return result;
    }

    public Map<Long, Integer> selectionVcfFieldPlotData(MgdbVcfFieldPlotRequest gvfpr, String token, boolean workWithSamples) throws Exception {
//      System.err.println(gvfpr.getVcfField() + " : " + gvfpr.getAllCallSetIds().stream().map(t -> t.size()).toList());

        long before = System.currentTimeMillis();

        String info[] = Helper.getInfoFromId(gvfpr.getVariantSetId(), 2);
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating plot data for " + gvfpr.getVcfField() +  " field regarding " + (gvfpr.getDisplayedVariantType() != null ? gvfpr.getDisplayedVariantType() + " " : "") + "variants on sequence " + gvfpr.getDisplayedSequence() });
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], AbstractTokenManager.readToken(gvfpr.getRequest()), false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());

        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gvfpr, true);
        Collection<BasicDBList> variantDataQueries = nTempVarCount == 0 ? varQueryWrapper.getVariantRunDataQueries() : varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();

        if (VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gvfpr, false).size() > 0 && nTempVarCount == 0) {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        final String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getNamespace().getCollectionName();
        final ConcurrentHashMap<Long, Integer> result = new ConcurrentHashMap<Long, Integer>();

        if (gvfpr.getDisplayedRangeMin() == null || gvfpr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gvfpr, nTempVarCount > 0 ? usedVarCollName : null))
                return result;

        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();
        final ProgressIndicator finalProgress = progress;

        HashSet<String> selectedMaterial = new HashSet<String>();
        List<List<String>> callsetIds = gvfpr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        String taskGroup = "vcfField_" + gvfpr.getVcfField() + "_" + progress.getProcessId();

        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(nTempVarCount > 0 ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gvfpr.getDisplayedRangeIntervalCount(), Arrays.asList(gvfpr.getDisplayedSequence()), gvfpr.getDisplayedVariantType(), gvfpr.getDisplayedRangeMin(), gvfpr.getDisplayedRangeMax(), nTempVarCount == 0 ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), 0);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (!finalProgress.isAborted())
                    {
                        List<String> variantsInInterval = mongoTemplate.getCollection(usedVarCollName).distinct("_id", intervalEntry.getValue(), String.class).into(new ArrayList<>());  // oddly, it is faster to run a pre-query than to use $lookup
                        final ArrayList<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();

                        BasicDBList matchList = new BasicDBList();
                        matchList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", variantsInInterval)));
                        if (nTempVarCount == 0 && !variantQueryDBList.isEmpty())
                            matchList.addAll(variantQueryDBList);
                        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", matchList)));

                        BasicDBObject projectStage = new BasicDBObject("$project", new BasicDBObject("r", new BasicDBObject("$let", new BasicDBObject()
                                .append("vars", new BasicDBObject("values", new BasicDBObject("$filter", new BasicDBObject()
                                        .append("input", individualOrSampleToCallSetListMap.values().stream().flatMap(List::stream).map(cs -> "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + cs.getId() + "." + SampleGenotype.SECTION_ADDITIONAL_INFO + "." + gvfpr.getVcfField()).toArray())
                                        .append("as", "val")
                                        .append("cond", new BasicDBObject("$ne", Arrays.asList("$$val", null)))
                                    )))
                                    .append("in", new BasicDBObject()
                                        .append("valSum", new BasicDBObject("$sum", "$$values"))
                                        .append("valCount", new BasicDBObject("$size", "$$values"))
                                    )
                                )));
                        pipeline.add(projectStage);

                        BasicDBObject groupStage = new BasicDBObject("$group", new BasicDBObject("_id", null)
                                .append("valSum", new BasicDBObject("$sum", new BasicDBObject("$toInt", "$r.valSum")))
                                .append("valCount", new BasicDBObject("$sum", "$r.valCount"))
                        );
                        pipeline.add(groupStage);

                        Iterator<Document> it = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).aggregate(pipeline).iterator();
                        if (!it.hasNext())
                            result.put(intervalEntry.getKey(), 0);
                        else {
                            Document resultDoc = it.next();
                            int totalCumulated = Double.valueOf(resultDoc.get("valSum").toString()).intValue();
                            result.put(intervalEntry.getKey(), totalCumulated == 0 ? 0 : (totalCumulated / Double.valueOf(resultDoc.get("valCount").toString()).intValue()));
                        }
                        finalProgress.setCurrentStepProgress((short) result.size() * 100 / gvfpr.getDisplayedRangeIntervalCount());
                        intervalEntry.setValue(null); // help GC
                    }
                }
            };

            threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(taskGroup, t)));
        }

        if (executor instanceof GroupedExecutor)
            ((GroupedExecutor) executor).shutdown(taskGroup);   // important to be sure that all tasks in the group are executed before the queue purges it
        else
            executor.shutdown();

        for (Future<Void> ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
            ttwf.get();

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionVcfFieldPlotData treated " + threadsToWaitFor.size() + " intervals on sequence " + gvfpr.getDisplayedSequence() + " between " + gvfpr.getDisplayedRangeMin() + " and " + gvfpr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<Long, Integer>(result);
    }

    public String igvData(MgdbChartRequest mdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(mdr.getVariantSetId());
        String processId = "igvViz_" + token;
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(processId))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(processId, new String[] { "Preparing data for IGV visualization" });
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();
        List<List<String>> allCallsetIDs = mdr.getAllCallSetIds();

        boolean fNoGenotypesRequested = allCallsetIDs.isEmpty() || (allCallsetIDs.size() == 1 && allCallsetIDs.get(0).isEmpty());
        Collection materialToExport = mdr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toList());
        List<Callset> callSetsToExport = workWithSamples
                ? ((Collection<ArrayList<Callset>>) MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, materialToExport).values()).stream().flatMap(Collection::stream).toList()
                : ((Collection<ArrayList<Callset>>) MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, materialToExport).values()).stream().flatMap(Collection::stream).toList();

        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        MongoCollection<Document> tempVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], token, false, false, false);
        boolean fWorkingOnTempColl = tempVarColl.countDocuments() > 0;

        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(mdr, true);
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (!fWorkingOnTempColl ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

        MongoCollection<Document> collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(fWorkingOnTempColl ? tempVarColl.getNamespace().getCollectionName() : mongoTemplate.getCollectionName(fNoGenotypesRequested ? VariantData.class : VariantRunData.class));
        StringBuffer sb = new StringBuffer();

        if (fNoGenotypesRequested) {    // simplest case where we're not returning genotypes: querying on variants collection will be faster
            Map<String, Integer> individualPositions = IExportHandler.buildIndividualPositions(callSetsToExport, workWithSamples);

            String header = "variant\talleles\tchrom\tpos";
            sb.append(header);
            for (String individual : individualPositions.keySet())
                sb.append("\t" + individual);
            sb.append("\n");

            MongoCursor<VariantData> varIt = collWithPojoCodec.find(variantQueryForTargetCollection, VariantData.class).projection(new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELES, 1).append(Assembly.getThreadBoundVariantRefPosPath(), 1)).iterator();
            while (varIt.hasNext()) {
                VariantData variant = varIt.next();
                ReferencePosition rp = variant.getReferencePosition(Assembly.getThreadBoundAssembly());
                sb.append(variant.getId() + "\t" + StringUtils.join(variant.getKnownAlleles(), "/") + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + (rp == null ? 0 : rp.getStartSite()) + "\n");
            }
        }
        else {
            final Map<Integer, String> callSetIdToBioEntityMap = new HashMap<>();
            for (Callset cs : callSetsToExport)
                callSetIdToBioEntityMap.put(cs.getId(), workWithSamples ? cs.getSampleId() : cs.getIndividual());

            Map<String, Collection<String>> bioEntitiesByPop = new HashMap<>();
            Map<String, HashMap<String, Float>> annotationFieldThresholdsByPop = new HashMap<>();
            List<List<String>> callsetIds = mdr.getAllCallSetIds();
            for (int i = 0; i < callsetIds.size(); i++) {
                bioEntitiesByPop.put(mdr.getGroupName(i), callsetIds.get(i).isEmpty() /* no selection means all selected */ ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));
                annotationFieldThresholdsByPop.put(mdr.getGroupName(i), mdr.getAnnotationFieldThresholds(i));
            }

            // count variants to display
            List<BasicDBObject> countPipeline = new ArrayList<>();
            if (!variantQueryForTargetCollection.isEmpty())
                countPipeline.add(new BasicDBObject("$match", variantQueryForTargetCollection));
            countPipeline.add(new BasicDBObject("$count", "count"));
            MongoCursor<Document> countCursor = (fWorkingOnTempColl ? collWithPojoCodec : mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class))).aggregate(countPipeline, Document.class).collation(IExportHandler.collationObj).iterator();

            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            HapMapExportHandler heh = (HapMapExportHandler) AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get("HAPMAP");
            heh.writeGenotypeFile(true, true, true, true, baos, info[0], mongoTemplate.findOne(new Query(Criteria.where("_id").is(Assembly.getThreadBoundAssembly())), Assembly.class), bioEntitiesByPop, workWithSamples, callSetIdToBioEntityMap, annotationFieldThresholdsByPop, progress, fWorkingOnTempColl ? tempVarColl.getNamespace().getCollectionName() : null, variantQueryForTargetCollection, countCursor.hasNext() ? ((Number) countCursor.next().get("count")).longValue() : 0, null, callSetsToExport);
            sb.append(baos.toString());
            progress.markAsComplete();
            LOG.debug("igvData processed range " + mdr.getDisplayedSequence() + ":" + mdr.getDisplayedRangeMin() + "-" + mdr.getDisplayedRangeMax() + " for " + new HashSet<>(callSetIdToBioEntityMap.values()).size() + " individuals in " + (System.currentTimeMillis() - before) / 1000f + "s");
        }

        return sb.toString();
    }
}