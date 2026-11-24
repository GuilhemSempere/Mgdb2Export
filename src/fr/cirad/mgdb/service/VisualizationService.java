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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
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
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.model.Collation;
import com.mongodb.client.model.CountOptions;
import com.mongodb.client.model.Projections;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.HapMapExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
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
import fr.cirad.tools.query.AdaptiveTaskSubmitter;
import fr.cirad.tools.query.GroupedExecutor;
import fr.cirad.tools.security.base.AbstractTokenManager;

/**
 * A service class responsible for generating dynamic visualization data
 * 
 * @author sempere
 */

@Component
public class VisualizationService {
    protected static final Logger LOG = Logger.getLogger(VisualizationService.class);

    /**
     * Number of concurrent chunk queries a chart starts from, and never goes below, when tuned (see createChartSubmitter()). Going lower
     * would mostly happen while other operations load the server, so that a chart would end up with too few queries in flight to get
     * its share of the server. It is kept low because charts on databases with many samples saturate the server early (each chunk
     * query then being costly), and starting above that point would bias the reference latency.
     */
    static final int MIN_CHART_CONCURRENCY = 10;

    /**
     * Creates the submitter through which a chart hands its interval queries to an executor obtained from MongoTemplateManager.getExecutor().
     * When maxQueryThreads_&lt;host&gt; is configured, all of them are queued at once, the shared pool's size bounding how many run at a time
     * and its threads being shared fairly between concurrent operations, as for searches (tuning that number proved slower there).
     * Otherwise, it is tuned from measured latency, starting from (and never going below) MIN_CHART_CONCURRENCY.
     */
    private static AdaptiveTaskSubmitter createChartSubmitter(ExecutorService executor, String taskGroup) {
        return executor instanceof GroupedExecutor ? AdaptiveTaskSubmitter.withoutLimit((GroupedExecutor) executor, taskGroup) : new AdaptiveTaskSubmitter(executor, taskGroup, MIN_CHART_CONCURRENCY, MIN_CHART_CONCURRENCY);
    }

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

        final ProgressIndicator finalProgress = progress;
        
        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());
        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), nTempVarCount == 0 ? variantQueryDBList : null, null);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), 0L);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (!finalProgress.isAborted()) {
                        long partialCount = mongoTemplate.getCollection(nTempVarCount == 0 ? mongoTemplate.getCollectionName(VariantData.class) : tmpVarColl.getNamespace().getCollectionName()).countDocuments(intervalEntry.getValue(), new CountOptions().collation(Helper.getPositionCollationIfFilteringOnSequence(intervalEntry.getValue())));
                        nTotalTreatedVariantCount.addAndGet((int) partialCount);
                        result.put(intervalEntry.getKey(), partialCount);
                        finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    }
                    else {
                        submitter.cancelAll();
                        LOG.debug("Cancelled query threads for process " + finalProgress.getProcessId());
                    }
                    intervalEntry.setValue(null); // help GC
                }
            };

            submitter.submit(t);
        }
        
        try {
            submitter.awaitCompletion();    // wait for all threads before moving to next phase
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

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
        BasicDBObject preCheckRangeQuery = new BasicDBObject();	// same criteria as the interval queries below, over the whole range
        if (sequences != null && !sequences.isEmpty())
            preCheckRangeQuery.put(refPosPathWithDot + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", sequences));
        if (variantType != null)
            preCheckRangeQuery.put(VariantData.FIELDNAME_TYPE, variantType);
        preCheckRangeQuery.put(refPosPathWithDot + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$gte", rangeMin).append("$lte", rangeMax));
        Collation precheckRangeCollation = Helper.getPositionCollationIfFilteringOnSequence(preCheckRangeQuery);
        long precheckRangeCount = precheckCollection == null ? 0 : (precheckRangeCollation == null ? precheckCollection.countDocuments(preCheckRangeQuery) : precheckCollection.countDocuments(preCheckRangeQuery, new CountOptions().collation(precheckRangeCollation)));
        ExecutorService preCheckExecutor = VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK <= 0 || precheckCollection == null || (float) precheckRangeCount / actualIntervalCount > VARIANTS_PER_INTERVAL_THRESHOLD_FOR_PRECHECK  ? null : Executors.newFixedThreadPool(Math.max(1, Runtime.getRuntime().availableProcessors() / 2));
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
                    boolean hasData = precheckCollection.find(query).collation(Helper.getPositionCollationIfFilteringOnSequence(query))
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

    /**
     * The variants a chart works on, and where it reads their genotypes: those of the temporary collection holding search results when
     * there is one, otherwise those of the main collection that match the search filters.
     *
     * Aggregations always run on variantRunData, without $lookup. On the main collection, each interval's query (sequence, position
     * range and search filters) is matched directly, with the position collation that lets it use the {sequence, position} index (which
     * was created with it: without it, the query would scan all variants of the position range, or the whole collection). Variants of a
     * temporary collection are first read from it (with that collation, for the same reason), then matched by ID without any collation:
     * an aggregation's collation applies to all its stages, and would prevent this match (or a $lookup on IDs) from using the variant ID
     * index, which has the simple collation.
     */
    private static final class ChartVariantSource {
        private final MongoCollection<Document> runDataColl;
        private final MongoCollection<Document> tempColl;	// null when working on the main collection
        private final BasicDBList variantQueryDBList;		// search filters, only applied to the main collection (a temporary one holds their results)
        private TreeMap<Long, List<Object>> tempVariantIdsByInterval;	// IDs of the temporary collection's variants, by interval start (see getChartIntervalQueries())

        private ChartVariantSource(MongoCollection<Document> runDataColl, MongoCollection<Document> tempColl, BasicDBList variantQueryDBList) {
            this.runDataColl = runDataColl;
            this.tempColl = tempColl;
            this.variantQueryDBList = variantQueryDBList;
        }

        /** @return the temporary collection's name, or null when working on the main collection */
        String temporaryCollectionName() {
            return tempColl == null ? null : tempColl.getNamespace().getCollectionName();
        }

        /** @return the IDs of the temporary collection's variants that fall in the interval starting at intervalStart, or null when working on the main collection */
        List<Object> temporaryVariantIds(Long intervalStart) {
            return tempColl == null ? null : tempVariantIdsByInterval.getOrDefault(intervalStart, Collections.emptyList());
        }

        /**
         * Runs stages on the run data of an interval's variants.
         *
         * @param intervalQuery the interval's query (see getIntervalQueries())
         * @param tempVariantIds the interval's variant IDs as returned by temporaryVariantIds()
         * @param scope the projects and runs involved, whose run data alone is matched
         * @param fBiallelicOnly whether to leave out variants with more than 2 known alleles
         * @param stages the stages to run on the matched variantRunData documents
         */
        AggregateIterable<Document> aggregate(BasicDBObject intervalQuery, List<Object> tempVariantIds, ChartRunScope scope, boolean fBiallelicOnly, List<BasicDBObject> stages) {
            BasicDBObject match = tempVariantIds == null ? new BasicDBObject(intervalQuery) : new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", tempVariantIds));
            if (scope.runFilter != null)
                match = new BasicDBObject("$and", Arrays.asList(match, scope.runFilter));
            if (fBiallelicOnly)
                match.append(VariantData.FIELDNAME_KNOWN_ALLELES + ".2", new BasicDBObject("$exists", false));

            List<BasicDBObject> pipeline = new ArrayList<>(stages.size() + 1);
            pipeline.add(new BasicDBObject("$match", match));
            pipeline.addAll(stages);
            AggregateIterable<Document> result = runDataColl.aggregate(pipeline).allowDiskUse(true);
            Collation collation = tempVariantIds == null ? Helper.getPositionCollationIfFilteringOnSequence(intervalQuery) : null;
            return collation == null ? result : result.collation(collation);
        }
    }

    /**
     * Opens a chart's variant source, or returns null (after setting the progress error) if the request filters on genotypes or annotations
     * while no search results are available.
     */
    private static ChartVariantSource openChartVariantSource(MgdbChartRequest gdr, String module, MongoTemplate mongoTemplate, ProgressIndicator progress) throws Exception {
        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gdr, true);
        Collection<BasicDBList> variantRunDataQueries = varQueryWrapper.getVariantRunDataQueries();

        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(module, AbstractTokenManager.readToken(gdr.getRequest()), false, false, false);
        boolean fGotTempVariants = tmpVarColl.find().projection(Projections.include("_id")).limit(1).first() != null;
        if (!fGotTempVariants && VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0) {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        return new ChartVariantSource(mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)), fGotTempVariants ? tmpVarColl : null,
                !fGotTempVariants && variantRunDataQueries.size() == 1 ? variantRunDataQueries.iterator().next() : null);
    }

    /**
     * Sets the request's displayed range to its default one if not specified.
     *
     * @return false (after setting the progress error) if no default range could be found
     */
    private boolean ensureChartRange(MgdbChartRequest gdr, ChartVariantSource source, ProgressIndicator progress, String sMethodName) throws Exception {
        if ((gdr.getDisplayedRangeMin() != null && gdr.getDisplayedRangeMax() != null) || findDefaultRangeMinMax(gdr, source.temporaryCollectionName()))
            return true;

        progress.setError(sMethodName + ": Unable to find default position range, either there are no results in the current selection for this sequence, or results are not in sync with interface filters (in which case, try re-applying filters).");
        return false;
    }

    /**
     * @return the chart's interval queries, keyed by interval start (see getIntervalQueries()), with null values for intervals known to hold
     * no variant. When working on a temporary collection, its variants are read at once, in a single pass over the displayed range, and
     * dispatched to the intervals by position (see ChartVariantSource.temporaryVariantIds()): resolving them interval by interval would double
     * the number of queries (and thus of round trips to the server, which dominate with thousands of small intervals on a remote server),
     * and that pass also tells which intervals are empty, without any precheck.
     */
    private HashMap<Long, BasicDBObject> getChartIntervalQueries(MgdbChartRequest gdr, ChartVariantSource source, MongoTemplate mongoTemplate) throws InterruptedException, ExecutionException {
        if (source.tempColl == null)
            return getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), source.variantQueryDBList, mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class)));

        HashMap<Long, BasicDBObject> intervalQueries = getIntervalQueries(gdr.getDisplayedRangeIntervalCount(), Arrays.asList(gdr.getDisplayedSequence()), gdr.getDisplayedVariantType(), gdr.getDisplayedRangeMin(), gdr.getDisplayedRangeMax(), null, null);
        TreeMap<Long, List<Object>> idsByInterval = new TreeMap<>();
        for (Long intervalStart : intervalQueries.keySet())
            idsByInterval.put(intervalStart, new ArrayList<>());

        String refPosPath = Assembly.getThreadBoundVariantRefPosPath(), startSitePath = refPosPath + "." + ReferencePosition.FIELDNAME_START_SITE;
        BasicDBObject rangeQuery = new BasicDBObject(refPosPath + "." + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", Arrays.asList(gdr.getDisplayedSequence())));
        if (gdr.getDisplayedVariantType() != null)
            rangeQuery.put(VariantData.FIELDNAME_TYPE, gdr.getDisplayedVariantType());
        rangeQuery.put(startSitePath, new BasicDBObject("$gte", gdr.getDisplayedRangeMin()).append("$lte", gdr.getDisplayedRangeMax()));
        FindIterable<Document> find = source.tempColl.find(rangeQuery).projection(Projections.include("_id", startSitePath));
        Collation collation = Helper.getPositionCollationIfFilteringOnSequence(rangeQuery);
        if (collation != null)
            find = find.collation(collation);
        for (Document variant : find) {
            Object position = variant;
            for (String pathPart : startSitePath.split("\\."))
                position = position instanceof Document ? ((Document) position).get(pathPart) : null;
            Long intervalStart = position instanceof Number ? idsByInterval.floorKey(((Number) position).longValue()) : null;	// intervals are contiguous: each one ends where the next one starts
            if (intervalStart != null)
                idsByInterval.get(intervalStart).add(variant.get("_id"));
        }

        for (Map.Entry<Long, List<Object>> intervalIds : idsByInterval.entrySet())
            if (intervalIds.getValue().isEmpty())
                intervalQueries.put(intervalIds.getKey(), null);
        source.tempVariantIdsByInterval = idsByInterval;
        return intervalQueries;
    }

    /** The individuals (or samples) a chart involves, by group, and their callsets in the involved projects */
    private static final class ChartMaterial {
        private final List<Collection<String>> groups = new ArrayList<>();
        private final List<Boolean> explicitlySelected = new ArrayList<>();	// false for a group with none selected, which then stands for all of the projects' ones
        private final TreeMap<String, List<Callset>> callsetsByEntity = new TreeMap<>();
        private int entityCount;	// number of distinct individuals (or samples), including any that has no callset in the projects

        /** @return a single group holding every involved individual (or sample) that has callsets in the projects */
        List<Collection<String>> pooled() {
            return Collections.singletonList(callsetsByEntity.keySet());
        }

        /**
         * @return the groups made disjoint, as required for per-group statistics: an individual (or sample) involved in several groups goes
         * to the last one that selected it explicitly, or if none did, to the last one involving it
         */
        List<Collection<String>> disjoint() {
            Map<String, Integer> groupByEntity = new HashMap<>();
            for (boolean fExplicitPass : new boolean[] {false, true})
                for (int i = 0; i < groups.size(); i++)
                    if (explicitlySelected.get(i) == fExplicitPass)
                        for (String entity : groups.get(i))
                            groupByEntity.put(entity, i);

            List<Collection<String>> result = new ArrayList<>();
            for (int i = 0; i < groups.size(); i++)
                result.add(new ArrayList<>());
            for (Map.Entry<String, Integer> entry : groupByEntity.entrySet())
                result.get(entry.getValue()).add(entry.getKey());
            return result;
        }
    }

    private static ChartMaterial getChartMaterial(MgdbChartRequest gdr, boolean workWithSamples) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(Integer::parseInt).toList();

        ChartMaterial material = new ChartMaterial();
        Collection<String> projectMaterial = null;	// fetched once, if needed
        HashSet<String> allEntities = new HashSet<>();
        for (List<String> groupCallSetIds : gdr.getAllCallSetIds()) {
            Collection<String> group;
            if (!groupCallSetIds.isEmpty())
                group = groupCallSetIds.stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toCollection(LinkedHashSet::new));
            else {
                if (projectMaterial == null)
                    projectMaterial = workWithSamples ? MgdbDao.getProjectSamples(info[0], projIDs) : MgdbDao.getProjectIndividuals(info[0], projIDs);
                group = projectMaterial;
            }
            material.groups.add(group);
            material.explicitlySelected.add(!groupCallSetIds.isEmpty());
            allEntities.addAll(group);
        }
        material.entityCount = allEntities.size();
        material.callsetsByEntity.putAll(workWithSamples ? MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, allEntities) : MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, allEntities));
        return material;
    }

    /** The projects and runs whose data a chart reads, i.e. those of the involved callsets */
    private static final class ChartRunScope {
        private BasicDBObject runFilter;	// restricts variantRunData documents to those runs, null if no callset is involved
        private boolean fMultipleRuns;		// whether a variant's genotypes may be spread over several documents
        private Set<Integer> allRunCallsetIds = new HashSet<>();	// every callset of these runs (whether involved or not)
    }

    private static ChartRunScope getChartRunScope(String module, MongoTemplate mongoTemplate, Map<String, List<Callset>> callsetsByEntity) throws Exception {
        ChartRunScope scope = new ChartRunScope();
        Map<Integer, Set<String>> involvedProjectRuns = callsetsByEntity.values().stream().flatMap(List::stream).collect(Collectors.groupingBy(Callset::getProjectId, Collectors.mapping(Callset::getRun, Collectors.toSet())));
        if (involvedProjectRuns.isEmpty())
            return scope;

        List<BasicDBObject> projectFilterList = new ArrayList<>();
        for (Map.Entry<Integer, Set<String>> projectRuns : involvedProjectRuns.entrySet()) {
            BasicDBObject projectFilter = new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, projectRuns.getKey());
            // if not all of this project's runs are involved, only match the required ones
            boolean fNotAllRunsNeeded = projectRuns.getValue().size() != mongoTemplate.findDistinct(new Query(Criteria.where("_id").is(projectRuns.getKey())), GenotypingProject.FIELDNAME_RUNS, GenotypingProject.class, String.class).size();
            projectFilterList.add(fNotAllRunsNeeded ? new BasicDBObject("$and", Arrays.asList(projectFilter, new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_RUNNAME, new BasicDBObject("$in", projectRuns.getValue())))) : projectFilter);
        }
        scope.runFilter = projectFilterList.size() == 1 ? projectFilterList.get(0) : new BasicDBObject("$or", projectFilterList);
        scope.fMultipleRuns = involvedProjectRuns.values().stream().mapToInt(Set::size).sum() > 1;

        for (Callset cs : MgdbDao.getCallsetsForProjectAndSamples(module, involvedProjectRuns.keySet(), null))
            if (involvedProjectRuns.containsKey(cs.getProjectId()) && involvedProjectRuns.get(cs.getProjectId()).contains(cs.getRun()))
                scope.allRunCallsetIds.add(cs.getId());
        return scope;
    }

    /** Field holding, for each group, the genotypes (strings) of its individuals (or samples) that have one (see buildGroupGenotypeStages()) */
    private static final String GROUP_GENOTYPES = "gts";

    /**
     * @param callsetIds IDs of callsets
     * @param fieldPath path of the wanted field within a callset's entry in sp (e.g. "gt")
     * @param scope the involved runs
     * @return an expression resolving to the value of that field for each of these callsets: when they are all the callsets of the involved
     * runs, it reads the whole sp object at once ($objectToArray), rather than looking up each callset's entry by name, which costs a scan of
     * its keys every time (i.e. grows with the square of the number of callsets)
     */
    private static Object callsetValuesExpr(Collection<Integer> callsetIds, String fieldPath, ChartRunScope scope) {
        if (!callsetIds.isEmpty() && callsetIds.size() == scope.allRunCallsetIds.size() && scope.allRunCallsetIds.containsAll(callsetIds))
            return new BasicDBObject("$map", new BasicDBObject("input", new BasicDBObject("$objectToArray", new BasicDBObject("$ifNull", Arrays.asList("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject()))))
                    .append("in", "$$this.v." + fieldPath));
        return callsetIds.stream().map(id -> "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + id + "." + fieldPath).collect(Collectors.toList());
    }

    /**
     * Builds the stages that turn variantRunData documents into one document per variant holding GROUP_GENOTYPES: for each group, the
     * genotypes of those of its individuals (or samples) that have one, missing data being left out, so that its size is the number of
     * called ones. Element order is irrelevant: every chart statistic is a sum or count over them.
     *
     * When several runs are involved, their documents are first merged (callset IDs being unique, so are the keys of their sp objects).
     * Individuals with a single callset contribute that callset's genotype (see callsetValuesExpr()). Those with several callsets contribute
     * the genotype their callsets agree on, none if they disagree.
     *
     * @param groups the individuals (or samples) of each group
     * @param callsetsByEntity callsets of each individual (or sample) in the involved projects
     */
    private static List<BasicDBObject> buildGroupGenotypeStages(ChartRunScope scope, List<Collection<String>> groups, Map<String, List<Callset>> callsetsByEntity) {
        List<BasicDBObject> stages = new ArrayList<>();
        if (scope.fMultipleRuns)
            stages.add(new BasicDBObject("$group", new BasicDBObject("_id", "$_id." + VariantRunDataId.FIELDNAME_VARIANT_ID)
                    .append(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$mergeObjects", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES))));

        String gtField = SampleGenotype.FIELDNAME_GENOTYPECODE;
        List<Object> groupGenotypeExprs = new ArrayList<>();
        for (Collection<String> group : groups) {
            List<Integer> singleCallsetIds = new ArrayList<>();
            List<Object> multipleCallsetGenotypes = new ArrayList<>();
            for (String entity : new TreeSet<>(group)) {
                List<Callset> callsets = callsetsByEntity.getOrDefault(entity, Collections.emptyList());
                if (callsets.size() == 1)
                    singleCallsetIds.add(callsets.get(0).getId());
                else if (callsets.size() > 1) {	// the genotype its callsets agree on, if any
                    Object calledGenotypes = new BasicDBObject("$filter", new BasicDBObject("input", callsetValuesExpr(callsets.stream().map(Callset::getId).toList(), gtField, scope))
                            .append("cond", new BasicDBObject("$ne", Arrays.asList("$$this", null))));
                    multipleCallsetGenotypes.add(letExpr(new BasicDBObject("obs", new BasicDBObject("$setUnion", Arrays.asList(calledGenotypes))),
                            condExpr(new BasicDBObject("$eq", Arrays.asList(new BasicDBObject("$size", "$$obs"), 1)), new BasicDBObject("$arrayElemAt", Arrays.asList("$$obs", 0)), null)));
                }
            }

            List<Object> parts = new ArrayList<>();
            if (!singleCallsetIds.isEmpty())
                parts.add(callsetValuesExpr(singleCallsetIds, gtField, scope));
            if (!multipleCallsetGenotypes.isEmpty())
                parts.add(multipleCallsetGenotypes);

            Object groupGenotypes = parts.isEmpty() ? new BasicDBObject("$literal", Collections.emptyList()) : (parts.size() == 1 ? parts.get(0) : new BasicDBObject("$concatArrays", parts));
            groupGenotypeExprs.add(new BasicDBObject("$filter", new BasicDBObject("input", groupGenotypes).append("cond", new BasicDBObject("$ne", Arrays.asList("$$this", null)))));
        }

        stages.add(new BasicDBObject("$project", new BasicDBObject("_id", 1).append(GROUP_GENOTYPES, groupGenotypeExprs)));
        return stages;
    }

    /* Statistics on a list of genotype strings (e.g. "0/1"), with shortcuts for the commonest diploid ones */

    /** @return the first (and, for pooled charts, only) group's genotypes */
    private static Object pooledGenotypes() {
        return new BasicDBObject("$arrayElemAt", Arrays.asList("$" + GROUP_GENOTYPES, 0));
    }

    private static Object genotypeAlleles(Object genotype) {
        return new BasicDBObject("$split", Arrays.asList(genotype, "/"));
    }

    private static Object diploidSwitch(Object genotype, Object homRef, Object het, Object homAlt, Object otherwise) {
        return new BasicDBObject("$switch", new BasicDBObject("branches", Arrays.asList(
                new BasicDBObject("case", new BasicDBObject("$eq", Arrays.asList(genotype, "0/0"))).append("then", homRef),
                new BasicDBObject("case", new BasicDBObject("$eq", Arrays.asList(genotype, "0/1"))).append("then", het),
                new BasicDBObject("case", new BasicDBObject("$eq", Arrays.asList(genotype, "1/1"))).append("then", homAlt)))
                .append("default", otherwise));
    }

    private static Object sumOverGenotypes(Object genotypes, Object perGenotypeExpr) {
        return new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", genotypes).append("as", "g").append("in", perGenotypeExpr)));
    }

    /** @return an expression summing the number of alternate alleles in genotypes (only meaningful for biallelic variants) */
    private static Object altAlleleCountExpr(Object genotypes) {
        Object altCount = new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject("input", genotypeAlleles("$$g")).append("cond", new BasicDBObject("$eq", Arrays.asList("$$this", "1")))));
        return sumOverGenotypes(genotypes, diploidSwitch("$$g", 0, 1, 2, altCount));
    }

    /** @return an expression summing the number of alleles (i.e. ploidy) in genotypes */
    private static Object alleleCountExpr(Object genotypes) {
        return sumOverGenotypes(genotypes, diploidSwitch("$$g", 2, 2, 2, new BasicDBObject("$size", genotypeAlleles("$$g"))));
    }

    /** @return an expression counting heterozygous genotypes */
    private static Object heterozygousCountExpr(Object genotypes) {
        Object fHet = new BasicDBObject("$gt", Arrays.asList(new BasicDBObject("$size", new BasicDBObject("$setUnion", Arrays.asList(genotypeAlleles("$$g")))), 1));
        return sumOverGenotypes(genotypes, condExpr(diploidSwitch("$$g", false, true, false, fHet), 1, 0));
    }

    private static Object letExpr(BasicDBObject vars, Object in) {
        return new BasicDBObject("$let", new BasicDBObject("vars", vars).append("in", in));
    }

    private static Object condExpr(Object ifExpr, Object thenExpr, Object elseExpr) {
        return new BasicDBObject("$cond", Arrays.asList(ifExpr, thenExpr, elseExpr));
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
        final ChartVariantSource source = openChartVariantSource(gdr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final ConcurrentHashMap<Long, Double> result = new ConcurrentHashMap<Long, Double>();
        if (!ensureChartRange(gdr, source, progress, "selectionFst"))
            return result;

        LOG.debug("Calculating Fst on " + info[0] + " with group sizes: " + gdr.getAllCallSetIds().stream().map(t -> t.size()).toList());
        ChartMaterial material = getChartMaterial(gdr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        final List<BasicDBObject> fstStages = buildFstStages(scope, material);
        final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gdr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    try {
                        List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                        Iterator<Document> it = tempVariantIds != null && tempVariantIds.isEmpty() ? Collections.emptyIterator() : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, true, fstStages).iterator();

                        if (progress.isAborted())
                            return;

                        /* Structure of a resulting document : {
                         *      _id: ...,
                         *      am: 1,  // highest allele index (biallelic variants only)
                         *      ps: [   // populations having data for this variant
                         *          {ss: 100, as: [ // sample size, alleles
                         *              {al: 0, af: 0.45, hf: 0.31},    // allele, allele frequency, heterozygote frequency
                         *              {al: 1, af: 0.55, hf: 0.31},
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
                        progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                        nTotalTreatedVariantCount.addAndGet(partialCount);
                    }
                    catch (MongoCommandException mce) {
                        progress.setError("Error during Fst calculation: " + mce.getErrorMessage());
                    }
                    intervalEntry.setValue(null); // help GC
                }
            };

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();    // wait for all threads before moving to next phase
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

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
        final ChartVariantSource source = openChartVariantSource(gdr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final ConcurrentHashMap<Long, Double> tajimaD = new ConcurrentHashMap<Long, Double>();
        final ConcurrentHashMap<Long, Double> segregatingSites = new ConcurrentHashMap<Long, Double>();
        if (!ensureChartRange(gdr, source, progress, "selectionTajimaD"))
            return new ArrayList<>();

        ChartMaterial material = getChartMaterial(gdr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        final List<BasicDBObject> tajimaDStages = buildTajimaDStages(gdr, mongoTemplate, scope, material);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gdr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
                segregatingSites.put(intervalEntry.getKey(), 0.0);
                tajimaD.put(intervalEntry.getKey(), Double.NaN);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                    Document chunk = tempVariantIds != null && tempVariantIds.isEmpty() ? null : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, true, tajimaDStages).first();
                    if (chunk != null) {
                        segregatingSites.put(intervalEntry.getKey(), ((Number) chunk.get(TJD_RES_SEGREGATINGSITES)).doubleValue());
                        tajimaD.put(intervalEntry.getKey(), chunk.getDouble(TJD_RES_TAJIMAD));
                    } else {
                        segregatingSites.put(intervalEntry.getKey(), 0.0);
                        tajimaD.put(intervalEntry.getKey(), Double.NaN);
                    }
                    progress.setCurrentStepProgress((short) segregatingSites.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();    // wait for all threads before moving to next phase
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionTajimaD treated " + submitter.getSubmittedCount() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");

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
        final ChartVariantSource source = openChartVariantSource(gdr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final Map<Long, Float> result = new ConcurrentHashMap<>();
        if (!ensureChartRange(gdr, source, progress, "selectionMaf"))
            return result;

        ChartMaterial material = getChartMaterial(gdr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        final List<BasicDBObject> mafStages = buildMafStages(scope, material);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());
        AtomicInteger totalVariantsProcessed = new AtomicInteger(0);

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gdr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                    Document chunk = tempVariantIds != null && tempVariantIds.isEmpty() ? null : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, true, mafStages).first();  // There's only one interval per query

                    if (chunk == null)
                        result.put(intervalEntry.getKey(), Float.NaN);
                    else {
                        int nVariantsInInterval = chunk.getInteger("n");
                        totalVariantsProcessed.addAndGet(nVariantsInInterval);
                        if (nVariantsInInterval == 0)
                            result.put(intervalEntry.getKey(), Float.NaN);
                        else {
                            float freq = ((Number) chunk.get("t")).floatValue() / nVariantsInInterval;
                            result.put(intervalEntry.getKey(), Math.min(freq, 100f - freq));
                        }
                    }

                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();    // wait for all threads before moving to next phase
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionMaf treated " + totalVariantsProcessed.get() + " variants on " + submitter.getSubmittedCount() + " intervals for sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return result;
    }

    private static final String FST_RES_SAMPLESIZE = "ss";
    private static final String FST_RES_ALLELEID = "al";
    private static final String FST_RES_ALLELEMAX = "am";
    private static final String FST_RES_ALLELEFREQUENCY = "af";
    private static final String FST_RES_HETEROZYGOUSFREQUENCY = "hf";
    private static final String FST_RES_ALLELES = "as";
    private static final String FST_RES_POPULATIONS = "ps";

    /** Fst needs per-allele frequencies (biallelic variants only) for each group, kept separate (one population each) */
    private static List<BasicDBObject> buildFstStages(ChartRunScope scope, ChartMaterial material) {
        List<BasicDBObject> stages = buildGroupGenotypeStages(scope, material.disjoint(), material.callsetsByEntity);

        // allele frequencies: 0 = ref, 1 = alt ; heterozygote frequency is the same for both alleles (biallelic)
        Object afRef = new BasicDBObject("$divide", Arrays.asList(new BasicDBObject("$subtract", Arrays.asList("$$totalCopies", "$$altCopies")), "$$totalCopies"));
        Object afAlt = new BasicDBObject("$divide", Arrays.asList("$$altCopies", "$$totalCopies"));
        Object hf = new BasicDBObject("$divide", Arrays.asList("$$heCount", "$$ss"));
        List<BasicDBObject> alleleEntries = Arrays.asList(
                new BasicDBObject(FST_RES_ALLELEID, 0).append(FST_RES_ALLELEFREQUENCY, afRef).append(FST_RES_HETEROZYGOUSFREQUENCY, hf),
                new BasicDBObject(FST_RES_ALLELEID, 1).append(FST_RES_ALLELEFREQUENCY, afAlt).append(FST_RES_HETEROZYGOUSFREQUENCY, hf));
        Object population = new BasicDBObject(FST_RES_SAMPLESIZE, "$$ss").append(FST_RES_ALLELES, alleleEntries);

        // populations with no data at all (ss == 0) are left out, so that selectionFst can tell whether at least 2 have some
        Object populationOrNull = letExpr(new BasicDBObject("ss", new BasicDBObject("$size", "$$gt")),
                condExpr(new BasicDBObject("$eq", Arrays.asList("$$ss", 0)), null,
                        letExpr(new BasicDBObject("altCopies", altAlleleCountExpr("$$gt")).append("totalCopies", alleleCountExpr("$$gt")).append("heCount", heterozygousCountExpr("$$gt")), population)));

        Object populations = new BasicDBObject("$filter", new BasicDBObject()
                .append("input", new BasicDBObject("$map", new BasicDBObject("input", "$" + GROUP_GENOTYPES).append("as", "gt").append("in", populationOrNull)))
                .append("cond", new BasicDBObject("$ne", Arrays.asList("$$this", null))));
        stages.add(new BasicDBObject("$project", new BasicDBObject("_id", 1).append(FST_RES_ALLELEMAX, new BasicDBObject("$literal", 1)).append(FST_RES_POPULATIONS, populations)));

        LOG.debug("Fst stages " + stages);
        return stages;
    }

    private static final String TJD_S23_ALLELEFREQUENCY = "k";
    private static final String TJD_S24_FREQUENCYSUM = "kc";

    private static final String TJD_RES_SEGREGATINGSITES = "sg";
    private static final String TJD_RES_TAJIMAD = "tjd";

    /** Tajima's D needs alt / total allele counts (biallelic variants only), all groups being pooled */
    private static List<BasicDBObject> buildTajimaDStages(MgdbChartRequest gdr, MongoTemplate mongoTemplate, ChartRunScope scope, ChartMaterial material) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gdr.getVariantSetId());
        Collection<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();

        // sampleSize = total allele count per locus assuming no missing data = ploidy × #individuals
        List<GenotypingProject> genotypingProjects = mongoTemplate.find(new Query(Criteria.where("_id").in(projIDs)), GenotypingProject.class);
        Integer[] ploidyLevels = genotypingProjects.stream().map(GenotypingProject::getPloidyLevel).distinct().toArray(Integer[]::new);
        if (ploidyLevels.length > 1)
            throw new Exception("Inconsistent ploidy levels among projects " + info[1] + " in database " + info[0]);
        final int sampleSize = ploidyLevels[0] * material.entityCount;

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

        List<BasicDBObject> stages = buildGroupGenotypeStages(scope, material.pooled(), material.callsetsByEntity);

        // Segregating ⟺ 0 < altAlleles < totalAlleles. For segregating sites, k = p(1-p); else k = null.
        Object segregatingCond = new BasicDBObject("$and", Arrays.asList(
                new BasicDBObject("$gt", Arrays.asList("$$altAlleles", 0)),
                new BasicDBObject("$lt", Arrays.asList("$$altAlleles", "$$totalAlleles"))));
        Object pq = letExpr(new BasicDBObject("freq", new BasicDBObject("$divide", Arrays.asList("$$altAlleles", "$$totalAlleles"))),
                new BasicDBObject("$multiply", Arrays.asList("$$freq", new BasicDBObject("$subtract", Arrays.asList(1, "$$freq")))));
        Object kExpr = letExpr(new BasicDBObject("gt", pooledGenotypes()),
                letExpr(new BasicDBObject("altAlleles", altAlleleCountExpr("$$gt")).append("totalAlleles", alleleCountExpr("$$gt")), condExpr(segregatingCond, pq, null)));
        stages.add(new BasicDBObject("$project", new BasicDBObject(TJD_S23_ALLELEFREQUENCY, kExpr)));

        // ── Group: sum p(1-p) over segregating sites, count segregating sites ──────────
        BasicDBObject outputGroup = new BasicDBObject();
        outputGroup.put("_id", 1);
        outputGroup.put(TJD_S24_FREQUENCYSUM, new BasicDBObject("$sum", new BasicDBObject("$ifNull", Arrays.asList("$" + TJD_S23_ALLELEFREQUENCY, 0))));
        outputGroup.put(TJD_RES_SEGREGATINGSITES, new BasicDBObject("$sum", new BasicDBObject("$cond", Arrays.asList(new BasicDBObject("$eq", Arrays.asList("$" + TJD_S23_ALLELEFREQUENCY, null)), 0, 1))));
        stages.add(new BasicDBObject("$group", outputGroup));

        // ── Compute Tajima's D ─────────────────────────────────────────────────────────
        // The $group above produces a document as soon as one variant is matched, even if none is segregating. With sg=0 the variance
        // term (e1*sg + e2*sg*(sg-1)) is 0, and MongoDB rejects a division by zero: report NaN instead (via $literal, since it can't be
        // produced by a computed $divide), consistent with the caller's chunk==null -> NaN handling for empty intervals.
        Object varianceTerm = new BasicDBObject("$sqrt", new BasicDBObject("$abs", new BasicDBObject("$add", Arrays.asList(
                new BasicDBObject("$multiply", Arrays.asList(e1, "$" + TJD_RES_SEGREGATINGSITES)),
                new BasicDBObject("$multiply", Arrays.asList(
                        new BasicDBObject("$multiply", Arrays.asList(e2, "$" + TJD_RES_SEGREGATINGSITES)),
                        new BasicDBObject("$subtract", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, 1))
                ))
        ))));

        Object tajimaDExpr = new BasicDBObject("$cond", new BasicDBObject()
                .append("if", new BasicDBObject("$eq", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, 0)))
                .append("then", new BasicDBObject("$literal", Double.NaN))
                .append("else", new BasicDBObject("$divide", Arrays.asList(
                        new BasicDBObject("$subtract", Arrays.asList(
                                new BasicDBObject("$divide", Arrays.asList(
                                        new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S24_FREQUENCYSUM, 2*sampleSize)),
                                        sampleSize - 1
                                )),
                                new BasicDBObject("$divide", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, a1))
                        )),
                        varianceTerm
                ))));

        BasicDBObject finalProject = new BasicDBObject();
        finalProject.put(TJD_RES_SEGREGATINGSITES, 1);
        finalProject.put(TJD_RES_TAJIMAD, tajimaDExpr);
        stages.add(new BasicDBObject("$project", finalProject));

        LOG.debug("Tajima's D stages " + stages);
        return stages;
    }

    /**
     * MAF (biallelic variants only) pools all groups into a single population. For each interval, sums the per-variant (folded) MAF as a
     * percentage ("t"), and counts the variants that have any genotype ("n").
     */
    private static List<BasicDBObject> buildMafStages(ChartRunScope scope, ChartMaterial material) {
        List<BasicDBObject> stages = buildGroupGenotypeStages(scope, material.pooled(), material.callsetsByEntity);

        // f = 100 × altAlleles / totalAlleles, null if no allele was called
        Object freqExpr = letExpr(new BasicDBObject("gt", pooledGenotypes()),
                letExpr(new BasicDBObject("altAlleles", altAlleleCountExpr("$$gt")).append("totalAlleles", alleleCountExpr("$$gt")),
                        condExpr(new BasicDBObject("$eq", Arrays.asList("$$totalAlleles", 0)), null,
                                new BasicDBObject("$divide", Arrays.asList(new BasicDBObject("$multiply", Arrays.asList("$$altAlleles", 100)), "$$totalAlleles")))));
        stages.add(new BasicDBObject("$project", new BasicDBObject("f", freqExpr)));

        Object fCalled = new BasicDBObject("$ne", Arrays.asList("$f", null));
        stages.add(new BasicDBObject("$group", new BasicDBObject("_id", null)
                .append("n", new BasicDBObject("$sum", condExpr(fCalled, 1, 0)))
                .append("t", new BasicDBObject("$sum", condExpr(fCalled, condExpr(new BasicDBObject("$lte", Arrays.asList("$f", 50)), "$f", new BasicDBObject("$subtract", Arrays.asList(100, "$f"))), 0)))));

        LOG.debug("MAF stages " + stages);
        return stages;
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
        final ChartVariantSource source = openChartVariantSource(gdr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final Map<Long, Float> result = new ConcurrentHashMap<>();
        if (!ensureChartRange(gdr, source, progress, "selectionMissingData"))
            return result;

        ChartMaterial material = getChartMaterial(gdr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        final List<BasicDBObject> missingDataStages = buildMissingDataStages(scope, material);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gdr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                    Document chunk = tempVariantIds != null && tempVariantIds.isEmpty() ? null : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, false, missingDataStages).first();  // There's only one interval per query
                    if (chunk == null)
                    	result.put(intervalEntry.getKey(), Float.NaN);
                    else {
	                    long totalNonMissing = ((Number) chunk.get("totalNonMissing")).longValue();
	                    long totalPossibleCalls = ((Number) chunk.get("totalPossibleCalls")).longValue();
	                    long totalMissingCalls  = totalPossibleCalls - totalNonMissing;
	                    result.put(intervalEntry.getKey(), totalPossibleCalls == 0 ? 0 : (float) totalMissingCalls / totalPossibleCalls * 100);
                    }
                    progress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    intervalEntry.setValue(null); // help GC
                }
            };

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionMissingData treated " + submitter.getSubmittedCount() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<>(result);
    }

    /**
     * Missing data (whatever the ploidy or number of alleles) pools all groups into a single population: sums, over the interval's variants,
     * the number of individuals (or samples) that have a genotype, and compares it with the number of possible calls, i.e. variants ×
     * individuals (including any that has no callset in the projects, which is always missing).
     */
    private static List<BasicDBObject> buildMissingDataStages(ChartRunScope scope, ChartMaterial material) {
        List<BasicDBObject> stages = buildGroupGenotypeStages(scope, material.pooled(), material.callsetsByEntity);
        stages.add(new BasicDBObject("$group", new BasicDBObject("_id", null)
                .append("totalNonMissing", new BasicDBObject("$sum", new BasicDBObject("$size", pooledGenotypes())))
                .append("variantCount", new BasicDBObject("$sum", 1))));
        stages.add(new BasicDBObject("$project", new BasicDBObject("totalNonMissing", 1)
                .append("totalPossibleCalls", new BasicDBObject("$multiply", Arrays.asList("$variantCount", material.entityCount)))));

        LOG.debug("Missing data stages: " + stages);
        return stages;
    }

    private static final String HET_S16_SAMPLESIZE = "ss";
    private static final String HET_RES_HETEROZYGOSITY = "he";
    private static final String HET_RES_POPULATIONS = "ps";

    /**
     * Heterozygosity (whatever the ploidy or number of alleles) pools all groups into a single population: for each variant, the proportion
     * of heterozygous genotypes in that population, as {ps: [{ss: <called genotypes>, he: <proportion>}]}, ps being empty if none was called.
     */
    private static List<BasicDBObject> buildHeterozygosityStages(ChartRunScope scope, ChartMaterial material) {
        List<BasicDBObject> stages = buildGroupGenotypeStages(scope, material.pooled(), material.callsetsByEntity);
        Object populations = letExpr(new BasicDBObject("gt", pooledGenotypes()),
                letExpr(new BasicDBObject("ss", new BasicDBObject("$size", "$$gt")),
                        condExpr(new BasicDBObject("$eq", Arrays.asList("$$ss", 0)), new BasicDBObject("$literal", Collections.emptyList()),
                                Arrays.asList(new BasicDBObject(HET_S16_SAMPLESIZE, "$$ss").append(HET_RES_HETEROZYGOSITY, new BasicDBObject("$divide", Arrays.asList(heterozygousCountExpr("$$gt"), "$$ss")))))));
        stages.add(new BasicDBObject("$project", new BasicDBObject("_id", 1).append(HET_RES_POPULATIONS, populations)));

        LOG.debug("Heterozygosity stages: " + stages);
        return stages;
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
        final ChartVariantSource source = openChartVariantSource(gdr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final Map<Long, Float> result = new ConcurrentHashMap<>();
        if (!ensureChartRange(gdr, source, progress, "selectionHeterozygosity"))
            return result;

        ChartMaterial material = getChartMaterial(gdr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        final List<BasicDBObject> heterozygosityStages = buildHeterozygosityStages(scope, material);

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, progress.getProcessId());

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gdr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), Float.NaN);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (progress.isAborted())
                        return;

                    List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                    Iterator<Document> it = tempVariantIds != null && tempVariantIds.isEmpty() ? Collections.emptyIterator() : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, false, heterozygosityStages).iterator();

                    if (progress.isAborted())
                        return;

                    // Collect all heterozygosity values across all variants and populations in this interval
                    List<Double> allHetValues = new ArrayList<>();
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

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, progress.getProcessId());
        }

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.info("selectionHeterozygosity treated " + submitter.getSubmittedCount() + " intervals on sequence " + gdr.getDisplayedSequence() + " between " + gdr.getDisplayedRangeMin() + " and " + gdr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");

        progress.markAsComplete();

        return new TreeMap<>(result);
    }

    public Map<Long, Integer> selectionVcfFieldPlotData(MgdbVcfFieldPlotRequest gvfpr, String token, boolean workWithSamples) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gvfpr.getVariantSetId());
        final ProgressIndicator progress = Optional.ofNullable(ProgressIndicator.get(token))
		    .orElseGet(() -> {
		    	ProgressIndicator pi = new ProgressIndicator(token, new String[] { "Calculating plot data for " + gvfpr.getVcfField() +  " field regarding " + (gvfpr.getDisplayedVariantType() != null ? gvfpr.getDisplayedVariantType() + " " : "") + "variants on sequence " + gvfpr.getDisplayedSequence() });
		        ProgressIndicator.registerProgressIndicator(pi);
		    	return pi;
		    });

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        final ChartVariantSource source = openChartVariantSource(gvfpr, info[0], mongoTemplate, progress);
        if (source == null)
            return null;

        final ConcurrentHashMap<Long, Integer> result = new ConcurrentHashMap<Long, Integer>();
        if (gvfpr.getDisplayedRangeMin() == null || gvfpr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gvfpr, source.temporaryCollectionName()))
                return result;

        // each involved callset's value of the field, in each run's document (summing over them needs no merging)
        ChartMaterial material = getChartMaterial(gvfpr, workWithSamples);
        final ChartRunScope scope = getChartRunScope(info[0], mongoTemplate, material.callsetsByEntity);
        List<Integer> callsetIds = material.callsetsByEntity.values().stream().flatMap(List::stream).map(Callset::getId).distinct().toList();
        Object values = callsetIds.isEmpty() ? new BasicDBObject("$literal", Collections.emptyList()) : callsetValuesExpr(callsetIds, SampleGenotype.SECTION_ADDITIONAL_INFO + "." + gvfpr.getVcfField(), scope);
        Object calledValues = new BasicDBObject("$filter", new BasicDBObject("input", values).append("cond", new BasicDBObject("$ne", Arrays.asList("$$this", null))));
        final List<BasicDBObject> stages = Arrays.asList(
                new BasicDBObject("$project", new BasicDBObject("r", letExpr(new BasicDBObject("values", calledValues),
                        new BasicDBObject("valSum", new BasicDBObject("$sum", "$$values")).append("valCount", new BasicDBObject("$size", "$$values"))))),
                new BasicDBObject("$group", new BasicDBObject("_id", null)
                        .append("valSum", new BasicDBObject("$sum", new BasicDBObject("$toInt", "$r.valSum")))
                        .append("valCount", new BasicDBObject("$sum", "$r.valCount"))));

        ExecutorService executor = MongoTemplateManager.getExecutor(info[0]);
        String taskGroup = "vcfField_" + gvfpr.getVcfField() + "_" + progress.getProcessId();
        final AdaptiveTaskSubmitter submitter = createChartSubmitter(executor, taskGroup);

        HashMap<Long, BasicDBObject> intervalQueries = getChartIntervalQueries(gvfpr, source, mongoTemplate);
        for (Map.Entry<Long, BasicDBObject> intervalEntry : intervalQueries.entrySet()) {
        	if (intervalEntry.getValue() == null) {
        		result.put(intervalEntry.getKey(), 0);
        		continue;
        	}

            Thread t = new Thread() {
                public void run() {
                    if (!progress.isAborted())
                    {
                        List<Object> tempVariantIds = source.temporaryVariantIds(intervalEntry.getKey());
                        Document resultDoc = tempVariantIds != null && tempVariantIds.isEmpty() ? null : source.aggregate(intervalEntry.getValue(), tempVariantIds, scope, false, stages).first();
                        if (resultDoc == null)
                            result.put(intervalEntry.getKey(), 0);
                        else {
                            int totalCumulated = Double.valueOf(resultDoc.get("valSum").toString()).intValue();
                            result.put(intervalEntry.getKey(), totalCumulated == 0 ? 0 : (totalCumulated / Double.valueOf(resultDoc.get("valCount").toString()).intValue()));
                        }
                        progress.setCurrentStepProgress((short) result.size() * 100 / gvfpr.getDisplayedRangeIntervalCount());
                        intervalEntry.setValue(null); // help GC
                    }
                }
            };

            submitter.submit(t);
        }

        try {
            submitter.awaitCompletion();    // wait for all threads before moving to next phase
        }
        finally {
            MongoTemplateManager.releaseExecutor(executor, taskGroup);
        }

        if (progress.isAborted())
            return null;

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionVcfFieldPlotData treated " + submitter.getSubmittedCount() + " intervals on sequence " + gvfpr.getDisplayedSequence() + " between " + gvfpr.getDisplayedRangeMin() + " and " + gvfpr.getDisplayedRangeMax() + " bp in " + (System.currentTimeMillis() - before)/1000f + "s");
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