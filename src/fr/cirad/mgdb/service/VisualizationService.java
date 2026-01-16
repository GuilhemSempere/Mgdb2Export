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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
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
import fr.cirad.model.MgdbDensityRequest;
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

    protected boolean findDefaultRangeMinMax(MgdbDensityRequest mdr, String tmpCollName /* if null, main variant coll is used*/) throws Exception
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
    
    public Map<Long, Long> selectionDensity(MgdbDensityRequest gdr, String processId) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());

        ProgressIndicator progress = new ProgressIndicator(processId, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "variant density on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionDensity: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filtersn).");
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
        LOG.debug("selectionDensity treated " + nTotalTreatedVariantCount.get() + " variants on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
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

    final int MIN_INTERVAL_SIZE = 10; // minimal interval size in bp

    public HashMap<Long, BasicDBObject> getIntervalQueries(
            int nIntervalCount,
            Collection<String> sequences,
            String variantType,
            long rangeMin,
            long rangeMax,
            BasicDBList variantQueryDBListToMerge,
            MongoCollection<Document> precheckCollection // can be null if we don't want to precheck
    ) throws InterruptedException, ExecutionException {

        HashMap<Long, BasicDBObject> result = new HashMap<>();
        if (rangeMax < rangeMin) return result;

        String refPosPathWithDot = Assembly.getThreadBoundVariantRefPosPath() + ".";

        long rangeSize = rangeMax - rangeMin;
        int actualIntervalCount = nIntervalCount;

        // Adjust interval count if requested interval would be too small
        double requestedIntervalSize = (double) rangeSize / Math.max(1, nIntervalCount);
        if (requestedIntervalSize < MIN_INTERVAL_SIZE && nIntervalCount > 1) {
            actualIntervalCount = (int) Math.max(1, Math.floor(rangeSize / MIN_INTERVAL_SIZE));
            LOG.debug(String.format("Adjusted interval count: %d -> %d (range: %d bp, requested interval size: %.1f bp < minimum %d bp)", nIntervalCount, actualIntervalCount, rangeSize, requestedIntervalSize, MIN_INTERVAL_SIZE));
        }
        actualIntervalCount = Math.max(1, actualIntervalCount);

        long intervalSize = (long) Math.ceil((double) rangeSize / actualIntervalCount);

        List<CompletableFuture<Void>> futures = new ArrayList<>();
        ExecutorService preCheckExecutor = 
        		precheckCollection == null || (float) precheckCollection.countDocuments(new BasicDBObject(refPosPathWithDot + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$gte", rangeMin).append("$lte", rangeMax))) / actualIntervalCount < 5 ? null
        			: Executors.newFixedThreadPool(Math.max(1, Runtime.getRuntime().availableProcessors() / 2));
        AtomicInteger keptQueryCount = new AtomicInteger(0);

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

    public Map<Long, Double> selectionFst(MgdbDensityRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());

        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Fst estimate on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionFst: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filtersn).");
                return result;
            }

        final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);
        final ProgressIndicator finalProgress = progress;

        List<BasicDBObject> baseQuery = buildFstQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();
        
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

            List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
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

//          System.err.println("submitting for " + progress.getProcessId());
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
    
    public List<Map<Long, Double>> selectionTajimaD(MgdbDensityRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());

        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Tajima's D on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionTajimaD: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filtersn).");
                return new ArrayList<>();
            }

        List<BasicDBObject> baseQuery = buildTajimaDQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
                segregatingSites.put(intervalEntry.getKey(), 0.0);
                tajimaD.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

        	List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
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
    
    public Map<Long, Float> selectionMaf(MgdbDensityRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        
        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "MAF on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionMaf: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filtersn).");
                return result;
            }

        List<BasicDBObject> baseQuery = buildMafQuery(gdr, useTempColl, workWithSamples);
        
        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            intervalEntry.getValue().append(VariantData.FIELDNAME_KNOWN_ALLELES + ".2", new BasicDBObject("$exists", false)); // exclude multi-allelic variants from MAF calculation
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));
            
            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    AggregateIterable<Document> queryResult = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(true);

                    Document chunk = queryResult.first();  // There's only one interval per query
                    
                    if (chunk == null)
                        result.put(intervalEntry.getKey(), Float.NaN);
                    else {
                        int nVariantsInInterval = chunk.getInteger("n");
                        if (nVariantsInInterval == 0)
                            result.put(intervalEntry.getKey(), Float.NaN);
                        else {
                            float freq = chunk.getDouble("t").floatValue() / nVariantsInInterval;
                            result.put(intervalEntry.getKey(), Math.min(freq, 100f - freq));
                        }
                    }
                    
                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
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
        LOG.debug("selectionMaf treated " + threadsToWaitFor.size() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
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

    private List<BasicDBObject> buildGenotypeDataQuery(MgdbDensityRequest gdr, boolean useTempColl, Map<String, List<Callset>> bioEntityToCallSetListMap, boolean keepPosition, boolean fGotMultiCallSetIndividuals, boolean workWithSamples) throws Exception {
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
            if (keepPosition)
                individualGroup.put(GENOTYPE_DATA_S7_POSITION, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S7_POSITION));
            pipeline.add(new BasicDBObject("$group", individualGroup));

            // Stage 11 : Weed out incoherent genotypes
            pipeline.add(new BasicDBObject("$match", new BasicDBObject(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$size", 1))));

            // Stage 12 : Group back by variant
            BasicDBObject variantGroup = new BasicDBObject();
            variantGroup.put("_id", "$_id");
            BasicDBObject spObject = new BasicDBObject();
            spObject.put("k", new BasicDBObject("$toString", "$_id." + GENOTYPE_DATA_S10_INDIVIDUALID));
            spObject.put("v", new BasicDBObject(TJD_S18_GENOTYPE, new BasicDBObject("$arrayElemAt", Arrays.asList("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, 0))));
            variantGroup.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$push", spObject));
            if (keepPosition)
                variantGroup.put(GENOTYPE_DATA_S7_POSITION, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S7_POSITION));
            pipeline.add(new BasicDBObject("$group", variantGroup));

            // Stage 13 : Convert back to sp object
            pipeline.add(new BasicDBObject("$project", new BasicDBObject("_id", "$_id." + FST_S22_VARIANTID).append(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$arrayToObject", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES))));
        }

        if (fWillNeedToMergeObjects || useTempColl) {  // merge genotypes into a single document per variant
            BasicDBObject projectMergeGroup = new BasicDBObject();
            projectMergeGroup.put("_id", "$_id" + (!useTempColl && !fGotMultiCallSetIndividuals ? "." + FST_S22_VARIANTID : ""));
            projectMergeGroup.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$mergeObjects", "$" + (useTempColl && !fGotMultiCallSetIndividuals ? GENOTYPE_DATA_S2_DATA + "." : "") + VariantRunData.FIELDNAME_SAMPLEGENOTYPES));
            pipeline.add(new BasicDBObject("$group", projectMergeGroup));
        }

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

    private List<BasicDBObject> buildFstQuery(MgdbDensityRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
//      System.err.println("Fst : " + gdr.getAllCallSetIds().stream().map(t -> t.size()).toList());
        
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        List<Collection<String>> selectedMaterial = new ArrayList<Collection<String>>();
        List<List<String>> ga4ghCallsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < ga4ghCallsetIds.size(); i++)
            selectedMaterial.add(ga4ghCallsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : ga4ghCallsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        for (Collection<String> group : selectedMaterial)
            if (!workWithSamples)
                individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, group));
            else 
                individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, group));

        boolean fGotMultiCallSetIndividuals = (individualOrSampleToCallSetListMap.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent());
        List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // Stage 14 : Get populations genotypes
        BasicDBList populationGenotypes = new BasicDBList();
        for (Collection<String> group : selectedMaterial)
            populationGenotypes.add(getFullPathToGenotypes(group, individualOrSampleToCallSetListMap, fGotMultiCallSetIndividuals));

        BasicDBObject projectGenotypes = new BasicDBObject(FST_S14_POPULATIONGENOTYPES, populationGenotypes);
        pipeline.add(new BasicDBObject("$project", projectGenotypes));

        // Stage 15 : Split by population
        BasicDBObject unwindPopulations = new BasicDBObject();
        unwindPopulations.put("path", "$" + FST_S14_POPULATIONGENOTYPES);
        unwindPopulations.put("includeArrayIndex", FST_S15_POPULATION);
        pipeline.add(new BasicDBObject("$unwind", unwindPopulations));

        // Stage 16 : Compute sample size
        BasicDBObject sampleSizeMapping = new BasicDBObject();
        sampleSizeMapping.put("input", "$" + FST_S14_POPULATIONGENOTYPES);
        sampleSizeMapping.put("in", new BasicDBObject("$cmp", Arrays.asList("$$this", null)));
        BasicDBObject addSampleSize = new BasicDBObject("$sum", new BasicDBObject("$map", sampleSizeMapping));
        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(FST_RES_SAMPLESIZE, addSampleSize)));

        // Stage 17 : Unwind by genotype
        pipeline.add(new BasicDBObject("$unwind", "$" + FST_S14_POPULATIONGENOTYPES));

        // Stage 18 : Eliminate missing genotypes
        BasicDBObject matchMissing = new BasicDBObject(FST_S14_POPULATIONGENOTYPES, new BasicDBObject("$ne", null));
        pipeline.add(new BasicDBObject("$match", matchMissing));

        // Stage 19 : Split genotype strings
        BasicDBObject projectSplitGenotypes = new BasicDBObject();
        projectSplitGenotypes.put(FST_S15_POPULATION, 1);
        projectSplitGenotypes.put(FST_RES_SAMPLESIZE, 1);
        BasicDBObject splitMapping = new BasicDBObject();
        splitMapping.put("input", new BasicDBObject("$split", Arrays.asList("$" + FST_S14_POPULATIONGENOTYPES, "/")));
        splitMapping.put("in", new BasicDBObject("$toInt", "$$this"));
        projectSplitGenotypes.put(FST_S19_GENOTYPE, new BasicDBObject("$map", splitMapping));
        pipeline.add(new BasicDBObject("$project", projectSplitGenotypes));

        // Stage 20 : Detect heterozygotes
        BasicDBList genotypeElements = new BasicDBList();  // TODO : Ploidy ?
        genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList("$" + FST_S19_GENOTYPE, 0)));
        genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList("$" + FST_S19_GENOTYPE, 1)));
        BasicDBObject addHeterozygote = new BasicDBObject(FST_S20_HETEROZYGOTE, new BasicDBObject("$ne", genotypeElements));
        pipeline.add(new BasicDBObject("$addFields", addHeterozygote));

        // Stage 21 : Unwind alleles
        pipeline.add(new BasicDBObject("$unwind", "$" + FST_S19_GENOTYPE));

        // Stage 22 : Group by allele
        BasicDBObject groupAllele = new BasicDBObject();
        BasicDBObject groupAlleleId = new BasicDBObject();
        groupAlleleId.put(FST_S22_VARIANTID, "$_id");
        groupAlleleId.put(FST_S22_POPULATIONID, "$" + FST_S15_POPULATION);
        groupAlleleId.put(FST_RES_ALLELEID, "$" + FST_S19_GENOTYPE);
        groupAllele.put("_id", groupAlleleId);
        groupAllele.put(FST_RES_SAMPLESIZE, new BasicDBObject("$first", "$" + FST_RES_SAMPLESIZE));
        groupAllele.put(FST_S22_ALLELECOUNT, new BasicDBObject("$sum", 1));
        groupAllele.put(FST_S22_HETEROZYGOUSCOUNT, new BasicDBObject("$sum", new BasicDBObject("$toInt", "$" + FST_S20_HETEROZYGOTE)));
        pipeline.add(new BasicDBObject("$group", groupAllele));

        // Stage 23 : Group by population
        BasicDBObject groupPopulation = new BasicDBObject();
        BasicDBObject groupPopulationId = new BasicDBObject();
        groupPopulationId.put(FST_S22_VARIANTID, "$_id." + FST_S22_VARIANTID);
        groupPopulationId.put(FST_S22_POPULATIONID, "$_id." + FST_S22_POPULATIONID);
        groupPopulation.put("_id", groupPopulationId);
        groupPopulation.put(FST_RES_SAMPLESIZE, new BasicDBObject("$first", "$" + FST_RES_SAMPLESIZE));
        groupPopulation.put(FST_RES_ALLELEMAX, new BasicDBObject("$max", "$_id." + FST_RES_ALLELEID));

        BasicDBObject groupPopulationAllele = new BasicDBObject();
        groupPopulationAllele.put(FST_RES_ALLELEID, "$_id." + FST_RES_ALLELEID);
        BasicDBObject alleleFrequencyOperation = new BasicDBObject("$divide", Arrays.asList("$" + FST_S22_ALLELECOUNT, new BasicDBObject("$multiply", Arrays.asList("$" + FST_RES_SAMPLESIZE, 2))));
        BasicDBObject hetFrequencyOperation = new BasicDBObject("$divide", Arrays.asList("$" + FST_S22_HETEROZYGOUSCOUNT, "$" + FST_RES_SAMPLESIZE));
        groupPopulationAllele.put(FST_RES_ALLELEFREQUENCY, alleleFrequencyOperation);
        groupPopulationAllele.put(FST_RES_HETEROZYGOUSFREQUENCY, hetFrequencyOperation);

        groupPopulation.put(FST_RES_ALLELES, new BasicDBObject("$push", groupPopulationAllele));
        pipeline.add(new BasicDBObject("$group", groupPopulation));

        // Stage 24 : Group by variant
        BasicDBObject groupVariant = new BasicDBObject();
        groupVariant.put("_id", "$_id." + FST_S22_VARIANTID);
        groupVariant.put(FST_RES_ALLELEMAX, new BasicDBObject("$max", "$" + FST_RES_ALLELEMAX));
        BasicDBObject groupVariantPopulation = new BasicDBObject();
        groupVariantPopulation.put(FST_RES_SAMPLESIZE, "$" + FST_RES_SAMPLESIZE);
        groupVariantPopulation.put(FST_RES_ALLELES, "$" + FST_RES_ALLELES);
        groupVariant.put(FST_RES_POPULATIONS, new BasicDBObject("$push", groupVariantPopulation));
        pipeline.add(new BasicDBObject("$group", groupVariant));

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

    private List<BasicDBObject> buildTajimaDQuery(MgdbDensityRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
//      System.err.println("Tajima : " + gdr.getAllCallSetIds().stream().map(t -> t.size()).toList());

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        HashSet<String> selectedMaterial = new HashSet<String>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else 
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));   

        final int sampleSize = 2*selectedMaterial.size();
        int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
        List<Long> intervalBoundaries = new ArrayList<Long>();
        for (int i = 0; i < gdr.getDisplayedRangeIntervalCount(); i++)
            intervalBoundaries.add(gdr.getDisplayedRangeMin() + (i*intervalSize));
        intervalBoundaries.add(gdr.getDisplayedRangeMax() + 1);

        double a1 = 0, a2 = 0;
        for (int i = 1; i < sampleSize; i++) {
            a1 += 1.0 / i;
            a2 += 1.0 / (i*i);
        }

        double b1 = (double)(sampleSize + 1) / (double)(3*(sampleSize - 1));
        double b2 = 2.0*(sampleSize*sampleSize + sampleSize + 3) / (9.0*sampleSize*(sampleSize - 1));
        double c1 = b1 - 1/a1;
        double c2 = b2 - (double)(sampleSize + 2)/(a1*sampleSize) + a2/(a1*a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1*a1 + a2);

        boolean fGotMultiCallSetIndividuals = (individualOrSampleToCallSetListMap.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent());
        List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualOrSampleToCallSetListMap, true, fGotMultiCallSetIndividuals, workWithSamples);
        String refPosPath = Assembly.getThreadBoundVariantRefPosPath();

        // Stage 14 : Get the genotypes needed
        DBObject gtArray;
        if (fGotMultiCallSetIndividuals)
        	gtArray = new BasicDBObject("$map", new BasicDBObject("input", new BasicDBObject("$objectToArray", "$" + GENOTYPE_DATA_S5_SPKEYVAL)).append("as", "cs").append("in", "$$cs.v." + TJD_S18_GENOTYPE));
        else
        	gtArray = getFullPathToGenotypes(selectedMaterial, individualOrSampleToCallSetListMap, fGotMultiCallSetIndividuals);
        BasicDBObject genotypeProjection = new BasicDBObject();
        genotypeProjection.put(refPosPath, 1);
        genotypeProjection.put(TJD_S14_GENOTYPES, gtArray);
        pipeline.add(new BasicDBObject("$project", genotypeProjection));

        // Stage 15 : Count non-null genotypes
        BasicDBObject sampleSizeMapping = new BasicDBObject();
        sampleSizeMapping.put("input", "$" + TJD_S14_GENOTYPES);
        sampleSizeMapping.put("in", new BasicDBObject("$cmp", Arrays.asList("$$this", null)));
        BasicDBObject addSampleSize = new BasicDBObject("$sum", new BasicDBObject("$map", sampleSizeMapping));
        pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(TJD_S15_SAMPLESIZE, addSampleSize)));

        // Stage 16 : Unwind individuals
        pipeline.add(new BasicDBObject("$unwind", "$" + TJD_S14_GENOTYPES));

        // Stage 17 : Eliminate missing genotypes
        BasicDBObject matchMissing = new BasicDBObject(TJD_S14_GENOTYPES, new BasicDBObject("$ne", null));
        pipeline.add(new BasicDBObject("$match", matchMissing));

        // Stage 18 : Split the genotype string
        BasicDBObject splitMapping = new BasicDBObject();
        splitMapping.put("input", new BasicDBObject("$split", Arrays.asList("$" + TJD_S14_GENOTYPES, "/")));
        splitMapping.put("in", new BasicDBObject("$toInt", "$$this"));
        BasicDBObject splitProjection = new BasicDBObject();
        splitProjection.put(refPosPath, 1);
        splitProjection.put(TJD_S18_GENOTYPE, new BasicDBObject("$map", splitMapping));
        splitProjection.put(TJD_S15_SAMPLESIZE, 1);
        pipeline.add(new BasicDBObject("$project", splitProjection));

        // Stage 19 : Unwind alleles
        pipeline.add(new BasicDBObject("$unwind", "$" + TJD_S18_GENOTYPE));

        // Stage 20 : Count the alleles
        BasicDBObject alleleGroupId = new BasicDBObject();
        alleleGroupId.put(TJD_S20_VARIANTID, "$_id");
        alleleGroupId.put(TJD_S20_ALLELEID, "$" + TJD_S18_GENOTYPE);
        BasicDBObject alleleGroup = new BasicDBObject();
        alleleGroup.put("_id", alleleGroupId);
        alleleGroup.put(VariantData.FIELDNAME_POSITIONS, new BasicDBObject("$first", "$" + refPosPath));
        alleleGroup.put(TJD_S20_ALLELECOUNT, new BasicDBObject("$sum", 1));
        alleleGroup.put(TJD_S15_SAMPLESIZE, new BasicDBObject("$first", "$" + TJD_S15_SAMPLESIZE));
        pipeline.add(new BasicDBObject("$group", alleleGroup));

        // Stage 21 : Group by variant, keeping only one of the two alleles
        BasicDBObject variantGroup = new BasicDBObject();
        variantGroup.put("_id", "$_id." + TJD_S20_VARIANTID);
        variantGroup.put(TJD_S20_ALLELECOUNT, new BasicDBObject("$first", "$" + TJD_S20_ALLELECOUNT));
        variantGroup.put(TJD_S21_NUMALLELES, new BasicDBObject("$sum", 1));
        variantGroup.put(TJD_S15_SAMPLESIZE, new BasicDBObject("$first", "$" + TJD_S15_SAMPLESIZE));
        variantGroup.put(VariantData.FIELDNAME_POSITIONS, new BasicDBObject("$first", "$" + VariantData.FIELDNAME_POSITIONS));
        pipeline.add(new BasicDBObject("$group", variantGroup));

        // Stage 22 : Keep only biallelic variants
        pipeline.add(new BasicDBObject("$match", new BasicDBObject(TJD_S21_NUMALLELES, 2)));

        // Stage 23 : Compute the average pairwise polymorphism for one variant
        BasicDBObject alleleFrequency = new BasicDBObject();
        alleleFrequency.put("vars", new BasicDBObject("freq",
                new BasicDBObject("$divide", Arrays.asList("$" + TJD_S20_ALLELECOUNT,
                    new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S15_SAMPLESIZE, 2))))));
        alleleFrequency.put("in", new BasicDBObject("$multiply", Arrays.asList("$$freq", new BasicDBObject("$subtract", Arrays.asList(1, "$$freq")))));  // p(1-p)
        BasicDBObject frequencyProjection = new BasicDBObject();
        frequencyProjection.put(VariantData.FIELDNAME_POSITIONS, 1);
        frequencyProjection.put(TJD_S23_ALLELEFREQUENCY, new BasicDBObject("$let", alleleFrequency));
        pipeline.add(new BasicDBObject("$project", frequencyProjection));

        // Stage 24 : Group by graph interval
        BasicDBObject outputGroup = new BasicDBObject();
        outputGroup.put("_id", 1);
        outputGroup.put(TJD_S24_FREQUENCYSUM, new BasicDBObject("$sum", "$" + TJD_S23_ALLELEFREQUENCY));
        outputGroup.put(TJD_RES_SEGREGATINGSITES, new BasicDBObject("$sum", 1));
        pipeline.add(new BasicDBObject("$group", outputGroup));

        // Stage 25 : Compute the Tajima's D value
        BasicDBObject finalProject = new BasicDBObject();
        finalProject.put(TJD_RES_SEGREGATINGSITES, 1);
        finalProject.put(TJD_RES_TAJIMAD, new BasicDBObject("$divide", Arrays.asList(
            new BasicDBObject("$subtract", Arrays.asList(
                new BasicDBObject("$divide", Arrays.asList(
                    new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S24_FREQUENCYSUM, 2*sampleSize)),
                    sampleSize - 1
                )),
                new BasicDBObject("$divide", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, a1))
            )),
            new BasicDBObject("$sqrt", new BasicDBObject("$abs", new BasicDBObject("$add", Arrays.asList(
                new BasicDBObject("$multiply", Arrays.asList(e1, "$" + TJD_RES_SEGREGATINGSITES)),
                new BasicDBObject("$multiply", Arrays.asList(
                    new BasicDBObject("$multiply", Arrays.asList(e2, "$" + TJD_RES_SEGREGATINGSITES)),
                    new BasicDBObject("$subtract", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, 1))
                ))
            ))))
        )));
        pipeline.add(new BasicDBObject("$project", finalProject));

        return pipeline;
    }

    private List<BasicDBObject> buildMafQuery(MgdbDensityRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
//      System.err.println("MAF : " + gdr.getAllCallSetIds().stream().map(t -> t.size()).toList());

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();
        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);

        List<GenotypingProject> genotypingProjects = mongoTemplate.find(new Query(Criteria.where("_id").in(projIDs)), GenotypingProject.class);
        Integer[] ploidyLevels = genotypingProjects.stream().map(pj -> pj.getPloidyLevel()).distinct().toArray(Integer[]::new);
        if (ploidyLevels.length > 1)
            throw new Exception("Inconsistent ploidy levels among projects " + info[1] + " in database " + info[0]);
        HashSet<String> selectedMaterial = new HashSet<String>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else 
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));   
        int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
        List<Long> intervalBoundaries = new ArrayList<Long>();
        for (int i = 0; i < gdr.getDisplayedRangeIntervalCount(); i++)
            intervalBoundaries.add(gdr.getDisplayedRangeMin() + (i*intervalSize));
        intervalBoundaries.add(gdr.getDisplayedRangeMax() + 1);

        boolean fGotMultiCallSetIndividuals = (individualOrSampleToCallSetListMap.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent());
        List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualOrSampleToCallSetListMap, true, fGotMultiCallSetIndividuals, workWithSamples);

        DBObject gtArray;
        if (fGotMultiCallSetIndividuals)
        	gtArray = new BasicDBObject("$map", new BasicDBObject("input", new BasicDBObject("$objectToArray", "$" + GENOTYPE_DATA_S5_SPKEYVAL)).append("as", "cs").append("in", "$$cs.v." + TJD_S18_GENOTYPE));
        else
        	gtArray = getFullPathToGenotypes(selectedMaterial, individualOrSampleToCallSetListMap, fGotMultiCallSetIndividuals);
        	
        BasicDBObject vars = new BasicDBObject(TJD_S18_GENOTYPE, gtArray);
        BasicDBObject in = new BasicDBObject();
        BasicDBObject subIn = new BasicDBObject();
        
        BasicDBObject inObj = new BasicDBObject("$add", Arrays.asList(1, new BasicDBObject("$cmp", Arrays.asList("$$g", ploidyLevels[0] == 1 ? "1" : "0/1"))));
        in.put("a", new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt").append("as", "g").append("in", inObj))));
        in.put("m", new BasicDBObject("$subtract", Arrays.asList(individualOrSampleToCallSetListMap.size(), new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt").append("as", "g").append("in", new BasicDBObject("$max", Arrays.asList(0, new BasicDBObject("$cmp", Arrays.asList("$$g", null))))))))));

        BasicDBList condList = new BasicDBList(), divideList = new BasicDBList();
        condList.add(new BasicDBObject("$eq", new Object[] {"$$m", individualOrSampleToCallSetListMap.size()}));
        condList.add(null);
        condList.add(new BasicDBObject("$subtract", new Object[] {individualOrSampleToCallSetListMap.size(), "$$m"}));
        divideList.add(new BasicDBObject("$multiply", new Object[] {"$$a", 50}));
        divideList.add(new BasicDBObject("$cond", condList));

        subIn.put("f", new BasicDBObject("$divide", divideList));
        subIn.put("n", new BasicDBObject("$cond", Arrays.asList(new BasicDBObject("$eq", new Object[] {"$$m", individualOrSampleToCallSetListMap.size()}), 0, 1)));

        BasicDBObject subVars = in;
        BasicDBObject subLet = new BasicDBObject("vars", subVars);
        subLet.put("in", subIn);
        in = new BasicDBObject("$let", subLet);

        BasicDBObject let = new BasicDBObject("vars", vars);
        let.put("in", in);
        
        Document projectionFields = new Document();
        projectionFields.put("r", new BasicDBObject("$let", let));
        pipeline.add(new BasicDBObject("$project", projectionFields));
        
        BasicDBObject intervalGroup = new BasicDBObject();  // prepare sum of frequencies + number of variants accounted for, so average can be calculated just by dividing
        intervalGroup.put("_id", null);
        intervalGroup.put("n", new BasicDBObject("$sum", "$r.n"));
        intervalGroup.put("t", new BasicDBObject("$sum", new BasicDBObject("$cond", Arrays.asList(new BasicDBObject("$lte", new Object[] {"$r.f", 50}), "$r.f", new BasicDBObject("$subtract", new Object[] {100, "$r.f"})))));
        pipeline.add(new BasicDBObject("$group", intervalGroup));

        return pipeline;
    }
    
    public Map<Long, Float> selectionMissingData(MgdbDensityRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        
        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "missing data percentage on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionMissingData: Unable to find default position range, either there are no results in the current selection for this sequence, or results are in sync with interface filters (in which case, try re-applying filtersn).");
                return result;
            }
        
        List<BasicDBObject> baseQuery = buildMissingDataQuery(gdr, useTempColl, workWithSamples);
        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        // for some reason only this type of query can be worth pre-filtering on lightwaeight collection to find empty intervals early
        MongoCollection<Document> precheckCollection = mongoTemplate.getCollection(useTempColl ? usedVarCollName : MongoTemplateManager.getMongoCollectionName(VariantData.class));
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, precheckCollection);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    Document chunk = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(true).first();  // There's only one interval per query
                    if (chunk == null)
                    	result.put(intervalEntry.getKey(), Float.NaN);
                    else {              
	                    int totalNonMissing = chunk.getInteger("totalNonMissing");
	                    int totalPossibleCalls = chunk.getInteger("totalPossibleCalls");
	                    int totalMissingCalls  = totalPossibleCalls - totalNonMissing;
	                    result.put(intervalEntry.getKey(), totalPossibleCalls == 0 ? 0 : (float) totalMissingCalls / totalPossibleCalls * 100);
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
    private List<BasicDBObject> buildMissingDataQuery(MgdbDensityRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
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
        List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualOrSampleToCallSetListMap, true,
                individualOrSampleToCallSetListMap.values().stream().anyMatch(spList -> spList.size() > 1),
                workWithSamples);

        // -------------------- $project: count non-missing individuals --------------------
        BasicDBObject filteredSpArray = new BasicDBObject("$filter",
                new BasicDBObject("input", new BasicDBObject("$objectToArray", "$sp"))
                        .append("as", "entry")
                        .append("cond", new BasicDBObject("$in", Arrays.asList("$$entry.k", allowedKeySet)))
        );

        // Count entries where gt != null
        BasicDBObject nonMissingCount = new BasicDBObject("$size",
                new BasicDBObject("$filter",
                        new BasicDBObject("input", filteredSpArray)
                                .append("as", "entry")
                                .append("cond", new BasicDBObject("$ne", Arrays.asList("$$entry.v.gt", null)))
                )
        );

        BasicDBObject projectStage = new BasicDBObject();
        projectStage.put("_id", 1);
        projectStage.put("nonMissingCount", nonMissingCount); // per-variant metric

        pipeline.add(new BasicDBObject("$project", projectStage));

        // -------------------- $group: aggregate across variants --------------------
        BasicDBObject groupStage = new BasicDBObject();
        groupStage.put("_id", null);
        groupStage.put("totalNonMissing", new BasicDBObject("$sum", "$nonMissingCount"));
        groupStage.put("variantCount", new BasicDBObject("$sum", 1));
        pipeline.add(new BasicDBObject("$group", groupStage));

        // -------------------- Final $project: compute totals and percentages --------------------
        BasicDBObject finalStage = new BasicDBObject();
        finalStage.put("totalPossibleCalls", new BasicDBObject("$multiply", Arrays.asList("$variantCount", totalIndividuals)));

        // total missing = totalIndividuals * variantCount − totalNonMissing
        finalStage.put("totalMissingCalls",
                new BasicDBObject("$subtract", Arrays.asList(
                        new BasicDBObject("$multiply", Arrays.asList("$variantCount", totalIndividuals)),
                        "$totalNonMissing"
                ))
        );

        finalStage.put("totalNonMissing", 1);
        pipeline.add(new BasicDBObject("$project", finalStage));

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
    private List<BasicDBObject> buildHeterozygosityQuery(MgdbDensityRequest gdr, boolean useTempColl, boolean workWithSamples) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        HashSet<String> selectedMaterial = new HashSet<String>();
        List<List<String>> callsetIds = gdr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++)
            selectedMaterial.addAll(callsetIds.get(i).isEmpty() ? (workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs)) : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<Callset>> individualOrSampleToCallSetListMap = new TreeMap<>();
        if (!workWithSamples)
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, selectedMaterial));
        else 
            individualOrSampleToCallSetListMap.putAll(MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, selectedMaterial));   

        boolean fGotMultiCallSetIndividuals = (individualOrSampleToCallSetListMap.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent());
        List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualOrSampleToCallSetListMap, false, fGotMultiCallSetIndividuals, workWithSamples);

        // Stage: Get populations genotypes        
        DBObject gtArray;
        if (fGotMultiCallSetIndividuals)
        	gtArray = new BasicDBObject("$map", new BasicDBObject("input", new BasicDBObject("$objectToArray", "$" + GENOTYPE_DATA_S5_SPKEYVAL)).append("as", "cs").append("in", "$$cs.v." + TJD_S18_GENOTYPE));
        else
        	gtArray = getFullPathToGenotypes(selectedMaterial, individualOrSampleToCallSetListMap, fGotMultiCallSetIndividuals);

        BasicDBObject projectGenotypes = new BasicDBObject();
        projectGenotypes.put("_id", 1);
        projectGenotypes.put(HET_S14_POPULATIONGENOTYPES, Arrays.asList(gtArray));	// this probably does not need to be an array of arrays
        pipeline.add(new BasicDBObject("$project", projectGenotypes));
        
        // Stage: Split by population
        BasicDBObject unwindPopulations = new BasicDBObject();
        unwindPopulations.put("path", "$" + HET_S14_POPULATIONGENOTYPES);
        unwindPopulations.put("includeArrayIndex", HET_S15_POPULATION);
        pipeline.add(new BasicDBObject("$unwind", unwindPopulations));

        // Stage: Compute sample size (non-null genotypes)
        BasicDBObject sampleSizeMapping = new BasicDBObject();
        sampleSizeMapping.put("input", "$" + HET_S14_POPULATIONGENOTYPES);
        sampleSizeMapping.put("in", new BasicDBObject("$cmp", Arrays.asList("$$this", null)));
        BasicDBObject addSampleSize = new BasicDBObject(HET_S16_SAMPLESIZE, new BasicDBObject("$sum", new BasicDBObject("$map", sampleSizeMapping)));
        pipeline.add(new BasicDBObject("$addFields", addSampleSize));

        // Stage: Unwind by genotype
        pipeline.add(new BasicDBObject("$unwind", "$" + HET_S14_POPULATIONGENOTYPES));

        // Stage: Eliminate missing genotypes
        BasicDBObject matchMissing = new BasicDBObject(HET_S14_POPULATIONGENOTYPES, new BasicDBObject("$ne", null));
        pipeline.add(new BasicDBObject("$match", matchMissing));

        // Stage: Check if heterozygote (for diploid: "0/1", for haploid: impossible so will be false)
        BasicDBObject addHeterozygote = new BasicDBObject();
        addHeterozygote.put(HET_S15_POPULATION, 1);
        addHeterozygote.put(HET_S16_SAMPLESIZE, 1);
        addHeterozygote.put(HET_S17_GENOTYPE, "$" + HET_S14_POPULATIONGENOTYPES);
        
        // Detect heterozygote: split genotype and check if alleles differ
        BasicDBObject splitGenotype = new BasicDBObject("$split", Arrays.asList("$" + HET_S14_POPULATIONGENOTYPES, "/"));
        BasicDBList genotypeElements = new BasicDBList();
        genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList(splitGenotype, 0)));
        genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList(splitGenotype, 1)));
        addHeterozygote.put(HET_S18_ISHETEROZYGOTE, new BasicDBObject("$ne", genotypeElements));
        
        pipeline.add(new BasicDBObject("$project", addHeterozygote));

        // Stage: Group by variant and population
        BasicDBObject groupVariantPopulation = new BasicDBObject();
        BasicDBObject groupVariantPopulationId = new BasicDBObject();
        groupVariantPopulationId.put("vi", "$_id");
        groupVariantPopulationId.put("pi", "$" + HET_S15_POPULATION);
        groupVariantPopulation.put("_id", groupVariantPopulationId);
        groupVariantPopulation.put(HET_S16_SAMPLESIZE, new BasicDBObject("$first", "$" + HET_S16_SAMPLESIZE));
        groupVariantPopulation.put(HET_RES_HETEROZYGOSITY, new BasicDBObject("$avg", new BasicDBObject("$toInt", "$" + HET_S18_ISHETEROZYGOTE)));
        pipeline.add(new BasicDBObject("$group", groupVariantPopulation));

        // Stage: Group by variant, collect populations
        BasicDBObject groupVariant = new BasicDBObject();
        groupVariant.put("_id", "$_id.vi");
        BasicDBObject populationData = new BasicDBObject();
        populationData.put(HET_S16_SAMPLESIZE, "$" + HET_S16_SAMPLESIZE);
        populationData.put(HET_RES_HETEROZYGOSITY, "$" + HET_RES_HETEROZYGOSITY);
        groupVariant.put(HET_RES_POPULATIONS, new BasicDBObject("$push", populationData));
        pipeline.add(new BasicDBObject("$group", groupVariant));

        return pipeline;
    }
    
    public Map<Long, Float> selectionHeterozygosity(MgdbDensityRequest gdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());

        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "heterozygosity on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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
                progress.setError("selectionHeterozygosity: Unable to find default position range, either there are no results in the current selection for this sequence, or results are not in sync with interface filters (in which case, try re-applying filtersn).");
                return result;
            }

        List<BasicDBObject> baseQuery = buildHeterozygosityQuery(gdr, useTempColl, workWithSamples);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();

        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), !useTempColl ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

        	List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
            windowQuery.set(0, new BasicDBObject("$match", intervalEntry.getValue()));

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    AggregateIterable<Document> queryResult = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(true);

                    if (progress.isAborted())
                        return;

                    // Collect all heterozygosity values across all variants and populations in this interval
                    List<Double> allHetValues = new ArrayList<>();
                    
                    Iterator<Document> it = queryResult.iterator();
                    while (it.hasNext()) {
                        Document variantResult = it.next();
                        List<Document> populations = variantResult.getList(HET_RES_POPULATIONS, Document.class);
                        
                        if (populations != null)
                            for (Document popData : populations) {
                                Object hetObj = popData.get(HET_RES_HETEROZYGOSITY);
                                if (hetObj != null) {
                                    double het = ((Number) hetObj).doubleValue();
                                    if (!Double.isNaN(het))
                                        allHetValues.add(het * 100);
                                }
                            }
                    }
                    
                    // Calculate average heterozygosity across all populations for this interval
                    result.put(intervalEntry.getKey(), allHetValues.isEmpty() ? Float.NaN : (float) allHetValues.stream().mapToDouble(Double::doubleValue).average().orElse(Double.NaN));
                    
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
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        ProgressIndicator progress = new ProgressIndicator(token, new String[] {"Calculating plot data for " + gvfpr.getVcfField() +  " field regarding " + (gvfpr.getDisplayedVariantType() != null ? gvfpr.getDisplayedVariantType() + " " : "") + "variants on sequence " + gvfpr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

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

        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gvfpr.getDisplayedRangeIntervalCount(), Arrays.asList(gvfpr.getDisplayedSequence()), gvfpr.getDisplayedVariantType(), gvfpr.getDisplayedRangeMin(), gvfpr.getDisplayedRangeMax(), nTempVarCount == 0 ? variantQueryDBList : null, null);
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

    public String igvData(MgdbDensityRequest mdr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(mdr.getVariantSetId());
        
        String processId = "igvViz_" + token;
        final ProgressIndicator progress = new ProgressIndicator(processId, new String[] {"Preparing data for visualization"});
        ProgressIndicator.registerProgressIndicator(progress);

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