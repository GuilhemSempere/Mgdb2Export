/*******************************************************************************
 * MGDB - Mongo Genotype DataBase
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
package fr.cirad.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import org.codehaus.jackson.annotate.JsonIgnoreProperties;
import org.ga4gh.methods.SearchVariantsRequest;

/**
 * ga4gh SearchVariantsRequest extended with extra option for filtering
 * @author petel, sempere, biggio
 */
@JsonIgnoreProperties(ignoreUnknown=true)
public class MgdbSearchVariantsRequest extends SearchVariantsRequest {

    private HttpServletRequest request;
    
    private String selectedVariantTypes = "";
    private String alleleCount = "";
    private String geneName = "";
    private String variantEffect = "";

    private List<String> groupName = new ArrayList<>();
    private List<String> gtPattern = new ArrayList<>();
    private List<HashMap<String, Float>> annotationFieldThresholds = new ArrayList<>();
    private List<Float> minMissingData = new ArrayList<>();
    private List<Float> maxMissingData = new ArrayList<>();
    private List<Float> minHeZ = new ArrayList<>();
    private List<Float> maxHeZ = new ArrayList<>();
    private List<Float> minMaf = new ArrayList<>();
    private List<Float> maxMaf = new ArrayList<>();
    private List<Integer> mostSameRatio = new ArrayList<>();
    private List<Integer> discriminate =  new ArrayList<>();
    private List<List<String>> additionalCallSetIds;    
    
    private String sortBy = "";
    private String sortDir = "asc";

    private int searchMode = 2;
    private boolean getGT = true;
    private boolean applyMatrixSizeLimit = true;
    
    private String selectedVariantIds = "";

	/** The Constant AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX. */
	static final public String AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX = "_ALL_"; // used to differentiate aggregation query with $and operator 

	/** The Constant AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX. */
	static final public String AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX = "_ATLO_";  // used to differentiate find query with $or operator

	/** The Constant AGGREGATION_QUERY_NEGATION_SUFFIX. */
	static final public String AGGREGATION_QUERY_NEGATION_SUFFIX = "_NEG_";    // used to indicate that the match operator should be negated in the aggregation query

	/** The Constant GENOTYPE_CODE_LABEL_ALL. */
	static final public String GENOTYPE_CODE_LABEL_ALL = "Any";

	/** The Constant GENOTYPE_CODE_LABEL_NOT_ALL_SAME. */
	static final public String GENOTYPE_CODE_LABEL_NOT_ALL_SAME = "Not all the same";

	/** The Constant GENOTYPE_CODE_LABEL_MOSTLY_SAME. */
	static final public String GENOTYPE_CODE_LABEL_MOSTLY_SAME = "All or mostly the same";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_DIFFERENT. */
	static final public String GENOTYPE_CODE_LABEL_ALL_DIFFERENT = "All different";

	/** The Constant GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT. */
	static final public String GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT = "Not all different";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF. */
	static final public String GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF = "All Homozygous Ref";

	/** The Constant GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF. */
	static final public String GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF = "Some Homozygous Ref";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR. */
	static final public String GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR = "All Homozygous Var";

	/** The Constant GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR. */
	static final public String GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR = "Some Homozygous Var";
    
    public MgdbSearchVariantsRequest() {
        super.setPageSize(100);
        super.setPageToken("0");
        super.setStart(-1L);
        super.setEnd(-1L);
    }

    public HttpServletRequest getRequest() {
        return request;
    }

    public void setRequest(HttpServletRequest request) {
        this.request = request;
    }
    
    public boolean shallApplyMatrixSizeLimit() {
        return applyMatrixSizeLimit;
    }

    public void setApplyMatrixSizeLimit(boolean applyMatrixSizeLimit) {
            this.applyMatrixSizeLimit = applyMatrixSizeLimit;
    }

    public int getSearchMode() {
        return searchMode;
    }

    public void setSearchMode(int searchMode) {
        this.searchMode = searchMode;
    }

    public String getSelectedVariantTypes() {
        return selectedVariantTypes;
    }

    public void setSelectedVariantTypes(String selectedVariantTypes) {
        this.selectedVariantTypes = selectedVariantTypes;
    }

    public Integer getNumberGroups() {
    	return Arrays.asList(gtPattern.size(), annotationFieldThresholds.size(), minMissingData.size(), maxMissingData.size(), minHeZ.size(), maxHeZ.size(), minMaf.size(), maxMaf.size(), mostSameRatio.size())
                .stream().mapToInt(Integer::intValue).max().orElse(0);
    }

    public List<List<String>> getAllCallSetIds() {
        return getNumberGroups() == 0 && getCallSetIds() == null || getCallSetIds().isEmpty() ? new ArrayList<>() : (additionalCallSetIds == null ? new ArrayList<>() {{ add(0, getCallSetIds()); }} : new ArrayList<>(additionalCallSetIds) {{ add(0, getCallSetIds()); }});
    }

    public boolean isGetGT() {
        return getGT;
    }

    public void setGetGT(boolean getGT) {
        this.getGT = getGT;
    }

    public String getGroupName(int nIndex) {
    	if (nIndex == 0 && groupName.isEmpty())
    		return null;	// no groups created, working with all individuals
		return groupName.get(nIndex);
	}

	public void setGroupNameWithIndex(String groupName, Integer index) {
		ensureSize(index);
		this.groupName.set(index, groupName);
	}

	public void setGroupName(List<String> groupName) {
		this.groupName = groupName;
	}

	public String getGtPattern(int nIndex) {
    	try {
    		return gtPattern.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return GENOTYPE_CODE_LABEL_ALL;
    	}
    }

    public void setGtPattern(List<String> gtPattern) {
        this.gtPattern = gtPattern;
    }
    
    public void setGtPatternWithIndex(String gtPattern, Integer index) {
    	ensureSize(index);
        this.gtPattern.set(index, gtPattern);
    }

    public void setDiscriminateWithIndex(Integer group, Integer index) {
    	ensureSize(index);
        this.discriminate.set(index, group);
    }

    private void ensureSize(Integer index) {
    	while (gtPattern.size() <= index) {
			gtPattern.add(GENOTYPE_CODE_LABEL_ALL);
			annotationFieldThresholds.add(new HashMap<>());
			minMissingData.add(0f);
			maxMissingData.add(100f);
			minHeZ.add(0f);
			maxHeZ.add(100f);
			minMaf.add(0f);
			maxMaf.add(50f);
			mostSameRatio.add(100);
			discriminate.add(null);
			if (gtPattern.size() < index)
				additionalCallSetIds.add(new ArrayList<>());
		}
	}
    
	public List<HashMap<String, Float>> getAnnotationFieldThresholds() {
   		return annotationFieldThresholds;
    }

	public HashMap<String, Float> getAnnotationFieldThresholds(int nIndex) {
    	try {
    		return annotationFieldThresholds.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return new HashMap<>();
    	}
    }

    public void setAnnotationFieldThresholds(List<HashMap<String, Float>> annotationFieldThresholds) {
        this.annotationFieldThresholds = annotationFieldThresholds;
    }

    public void setAnnotationFieldThresholdsWithIndex(HashMap<String, Float> annotationFieldThresholds, Integer index) {
    	ensureSize(index);
        this.annotationFieldThresholds.set(index, annotationFieldThresholds);
    }

    public Float getMinMissingData(int nIndex) {
    	try {
    		return minMissingData.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 0f;
    	}
    }

    public void setMinMissingData(List<Float> minMissingData) {
        this.minMissingData = minMissingData;
    }

    public void setMinMissingDataWithIndex(Float minMissingData, Integer index) {
    	ensureSize(index);
        this.minMissingData.set(index, minMissingData);
    }

    public Float getMaxMissingData(int nIndex) {
    	try {
    		return maxMissingData.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 100f;
    	}
    }

    public void setMaxMissingData(List<Float> maxMissingData) {
        this.maxMissingData = maxMissingData;
    }

    public void setMaxMissingDataWithIndex(Float maxMissingData, Integer index) {
    	ensureSize(index);
        this.maxMissingData.set(index, maxMissingData);
    }

    public Float getMinHeZ(int nIndex) {
    	try {
    		return minHeZ.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 0f;
    	}
    }

    public void setMinHeZ(List<Float> minHeZ) {
        this.minHeZ = minHeZ;
    }

    public void setMinHeZWithIndex(Float minHeZ, Integer index) {
    	ensureSize(index);
        this.minHeZ.set(index, minHeZ);
    }

    public Float getMaxHeZ(int nIndex) {
    	try {
    		return maxHeZ.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 100f;
    	}
    }

    public void setMaxHeZ(List<Float> maxHeZ) {
        this.maxHeZ = maxHeZ;
    }

    public void setMaxHeZWithIndex(Float maxHeZ, Integer index) {
    	ensureSize(index);
        this.maxHeZ.set(index, maxHeZ);
    }

    public Float getMinMaf(int nIndex) {
    	try {
    		return minMaf.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 0f;
    	}
    }

    public void setMinMaf(List<Float> minMaf) {
        this.minMaf = minMaf;
    }

    public void setMinMafWithIndex(Float minMaf, Integer index) {
    	ensureSize(index);
        this.minMaf.set(index, minMaf);
    }

    public Float getMaxMaf(int nIndex) {
    	try {
    		return maxMaf.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 50f;
    	}
    }

    public void setMaxMaf(List<Float> maxMaf) {
        this.maxMaf = maxMaf;
    }

    public void setMaxMafWithIndex(Float maxMaf, Integer index) {
    	ensureSize(index);
        this.maxMaf.set(index, maxMaf);
    }

    public String getAlleleCount() {
        return alleleCount;
    }

    public void setAlleleCount(String alleleCount) {
        this.alleleCount = alleleCount;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getVariantEffect() {
        return variantEffect;
    }

    public void setVariantEffect(String variantEffect) {
        this.variantEffect = variantEffect;
    }
    
    public int getMostSameRatio(int nIndex) {
    	try {
    		return mostSameRatio.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return 100;
    	}
    }

    public void setMostSameRatio(List<Integer> mostSameRatio) {
        this.mostSameRatio = mostSameRatio;
    }

    public void setMostSameRatioWithIndex(Integer mostSameRatio, Integer index) {
    	ensureSize(index);
        this.mostSameRatio.set(index, mostSameRatio);
    }

    public List<String> getAdditionalCallSetIds(int nIndex) {
    	try {
    		return additionalCallSetIds.get(nIndex);
    	}
    	catch (IndexOutOfBoundsException iobe) {
    		return new ArrayList<>();
    	}
    }

    public void setAdditionalCallSetIds(List<List<String>> additionalCallSetIds) {
        this.additionalCallSetIds = additionalCallSetIds;
    }

    public boolean isDiscriminate(int groupNumber) {
        return discriminate.get(groupNumber) != null || discriminate.contains(groupNumber + 1);
    }

	public List<Integer> getDiscriminate() {
		return discriminate;
	}

    public void setDiscriminate(List<Integer> discriminate) {
        this.discriminate = discriminate;
    }

    public String getSortBy() {
        return sortBy;
    }

    public void setSortBy(String sortBy) {
        this.sortBy = sortBy;
    }

    public String getSortDir() {
        return sortDir;
    }

    public void setSortDir(String sortDir) {
        this.sortDir = sortDir;
    }

    public String getSelectedVariantIds() {
        return selectedVariantIds;
    }

    public void setSelectedVariantIds(String selectedVariantIds) {
        this.selectedVariantIds = selectedVariantIds;
    }
}