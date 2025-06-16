/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
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
package fr.cirad.tools.mgdb;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;

import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;

import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.Run;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.MgdbDensityRequest;
import fr.cirad.model.MgdbSearchVariantsRequest;
import fr.cirad.tools.Helper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class VariantQueryBuilder.
 */
public class VariantQueryBuilder
{
    public static final Integer QUERY_IDS_CHUNK_SIZE = 100000;

    static public VariantQueryWrapper buildVariantDataQuery(MgdbSearchVariantsRequest gsvr, /*List<String> externallySelectedSeqs, */boolean fForBrowsing) throws Exception {
        String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(gsvr.getVariantSetId());
        Integer[] projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toArray(Integer[]::new);

        String actualSequenceSelection = gsvr.getReferenceName();
//        if (actualSequenceSelection == null || actualSequenceSelection.length() == 0) {
//            if (externallySelectedSeqs != null) {
//                actualSequenceSelection = StringUtils.join(externallySelectedSeqs, ";");
//            }
//        }
        List<String> selectedVariantTypes = gsvr.getSelectedVariantTypes().length() == 0 ? null : Arrays.asList(gsvr.getSelectedVariantTypes().split(";"));
        List<String> selectedSequences = Arrays.asList(actualSequenceSelection == null || actualSequenceSelection.length() == 0 ? new String[0] : actualSequenceSelection.split(";"));
        List<String> alleleCountList = gsvr.getAlleleCount().length() == 0 ? null : Arrays.asList(gsvr.getAlleleCount().split(";"));
        List<String> selectedVariantIds = gsvr.getSelectedVariantIds().length() == 0 ? null : Arrays.asList(gsvr.getSelectedVariantIds().split(";"));

        Collection<BasicDBList> queries = new ArrayList<>();
        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        BasicDBObject projectQueryForVar = null;
		BasicDBObject projectQueryForVrd = null;

        if (selectedVariantIds == null || fForBrowsing) {
            BasicDBList variantFeatureFilterList = new BasicDBList();
            
            // Filter on project
            if (Helper.estimDocCount(mongoTemplate, GenotypingProject.class) != 1) {
            	projectQueryForVar = new BasicDBObject(VariantData.FIELDNAME_RUNS + "." + Run.FIELDNAME_PROJECT_ID, new BasicDBObject("$in", projIDs));
            	projectQueryForVrd = new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, new BasicDBObject("$in", projIDs));
            }
            
            /* Step to match selected variant types */
            if (selectedVariantTypes != null && selectedVariantTypes.size() > 0) {
                BasicDBList orList1 = new BasicDBList();
                BasicDBObject orSelectedVariantTypesList = new BasicDBObject();
                for (String aSelectedVariantTypes : selectedVariantTypes) {
                    BasicDBObject orClause1 = new BasicDBObject(VariantData.FIELDNAME_TYPE, aSelectedVariantTypes);
                    orList1.add(orClause1);
                    orSelectedVariantTypesList.put("$or", orList1);
                }
                variantFeatureFilterList.add(orSelectedVariantTypesList);
            }
            
            String refPosPath = Assembly.getThreadBoundVariantRefPosPath();
           
            /* Step to match variants by position */
            LinkedHashSet<BasicDBObject> posOrSet = new LinkedHashSet<>();
            ArrayList<String> leftBounds = new ArrayList<>() {{ add(refPosPath + "." + ReferencePosition.FIELDNAME_START_SITE); }};
            Collection<String> variantTypes = MgdbDao.getVariantTypes(MongoTemplateManager.get(info[0]), projIDs);
            if (variantTypes.size() != 1 || !Type.SNP.toString().equals(variantTypes.iterator().next()))
            	leftBounds.add(refPosPath + "." + ReferencePosition.FIELDNAME_END_SITE);	// only add this part if the project contains non-SNP variants (otherwise it unnecessarily complexifies & slows down query execution)
            
            for (String leftBound : leftBounds) {
            	ArrayList<BasicDBObject> posAndList = new ArrayList<>();

                /* match variants that have a position included in the specified range */
                if (gsvr.getStart() != null && gsvr.getStart() != -1)
                	posAndList.add(new BasicDBObject(leftBound, new BasicDBObject("$gte", gsvr.getStart())));
                if (gsvr.getEnd() != null && gsvr.getEnd() != -1)
                	posAndList.add(new BasicDBObject(refPosPath + "." + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$lte", gsvr.getEnd())));
                
                /* match selected chromosomes */
                if (selectedSequences != null && selectedSequences.size() > 0)
                	posAndList.add(new BasicDBObject(refPosPath + "." + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", selectedSequences)));
                
                if (!posAndList.isEmpty())
                	posOrSet.add(new BasicDBObject("$and", posAndList));
            }
            if (!posOrSet.isEmpty())
            	variantFeatureFilterList.add(new BasicDBObject("$or", posOrSet));
            
            /* Step to match variants position range for visualization (differs from the above, which is for defining the subset of data Gigwa is currently working with: it they are contradictory it still makes sense and means user is trying to view variants outside the range selected in Gigwa) */
            if (MgdbDensityRequest.class.isAssignableFrom(gsvr.getClass())) {
                variantFeatureFilterList.add(new BasicDBObject(refPosPath + "." + ReferencePosition.FIELDNAME_SEQUENCE, ((MgdbDensityRequest) gsvr).getDisplayedSequence()));

                BasicDBObject posCrit = new BasicDBObject();
                Long min = ((MgdbDensityRequest) gsvr).getDisplayedRangeMin(), max = ((MgdbDensityRequest) gsvr).getDisplayedRangeMax();
                if (min != null && min != -1)
                    posCrit.put("$gte", min);
                if (max != null && max != -1)
                    posCrit.put("$lte", max);
                if (!posCrit.isEmpty())
                variantFeatureFilterList.add(new BasicDBObject(refPosPath + "." + ReferencePosition.FIELDNAME_START_SITE, posCrit));
            }
            
            /* Step to match selected number of alleles */
            if (alleleCountList != null) {
                BasicDBList orList3 = new BasicDBList();
                BasicDBObject orSelectedNumberOfAllelesList = new BasicDBObject();
                for (String aSelectedNumberOfAlleles : alleleCountList) {
                    int alleleNumber = Integer.parseInt(aSelectedNumberOfAlleles);
                    orList3.add(new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELES, new BasicDBObject("$size", alleleNumber)));
                    orSelectedNumberOfAllelesList.put("$or", orList3);
                }
                variantFeatureFilterList.add(orSelectedNumberOfAllelesList);
            }
            if (!variantFeatureFilterList.isEmpty())
                queries.add(variantFeatureFilterList);
        }
        else {    // filtering on variant IDs: we might need to split the query in order to avoid reaching a 16Mb document size
            int step = selectedVariantIds.size() / QUERY_IDS_CHUNK_SIZE;
            int r = selectedVariantIds.size() % QUERY_IDS_CHUNK_SIZE;
            if (r != 0)
                step++;

            for (int i = 0; i<step; i++) {
                List<String> subList = selectedVariantIds.subList(i*QUERY_IDS_CHUNK_SIZE, Math.min((i+1)*QUERY_IDS_CHUNK_SIZE, selectedVariantIds.size()));
                BasicDBList variantFeatureFilterList = new BasicDBList();
                variantFeatureFilterList.add(new BasicDBObject("_id", new BasicDBObject("$in", subList)));
                queries.add(variantFeatureFilterList);
            }
        }
        return new VariantQueryWrapper(queries, projectQueryForVar, projectQueryForVrd);
    }
    
    static public List<Integer> getGroupsForWhichToFilterOnGenotypingData(MgdbSearchVariantsRequest gsvr, boolean fConsiderFieldThresholds)
    {
        List<Integer> result = new ArrayList<>();
        for (int i = 0; i < gsvr.getNumberGroups(); i++)
            if (gsvr.isDiscriminate(i) || !gsvr.getGtPattern(i).equals(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL) || gsvr.getMinHeZ(i) > 0 || gsvr.getMaxHeZ(i) < 100 || gsvr.getMinMissingData(i) > 0 || gsvr.getMaxMissingData(i) < 100 || gsvr.getMinMaf(i) > 0 || gsvr.getMaxMaf(i) < 50)
                result.add(i);
        return result;
    }
    
    static public List<Integer> getGroupsForWhichToFilterOnGenotypingOrAnnotationData(MgdbSearchVariantsRequest gsvr, boolean fConsiderFielThresholds)
    {
        List<Integer> result = getGroupsForWhichToFilterOnGenotypingData(gsvr, fConsiderFielThresholds);
 
        if (result.size() == 0 && (gsvr.getGeneName().length() > 0 || gsvr.getVariantEffect().length() > 0))
            result.add(0);    // needed at least for filtering on annotation data or distinguish records according to project id

        /*FIXME: this should also force filtering on VRD in cases where only some runs of the selected project are involved*/
        return result;
    }
}