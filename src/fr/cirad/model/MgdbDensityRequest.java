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

/**
 *
 * @author petel, sempere
 */
public class MgdbDensityRequest extends MgdbSearchVariantsRequest {

    private String displayedSequence;
    private Long displayedRangeMin;
    private Long displayedRangeMax;
    private Integer displayedRangeIntervalCount;
    private String displayedVariantType;

	public MgdbDensityRequest(){
        super();
    }
    
    public String getDisplayedSequence() {
        return displayedSequence;
    }

    public void setDisplayedSequence(String displayedSequence) {
        this.displayedSequence = displayedSequence;
    }

    public Long getDisplayedRangeMin() {
        return displayedRangeMin;
    }

    public void setDisplayedRangeMin(Long displayedRangeMin) {
        this.displayedRangeMin = displayedRangeMin;
    }

    public Long getDisplayedRangeMax() {
        return displayedRangeMax;
    }

    public void setDisplayedRangeMax(Long displayedRangeMax) {
        this.displayedRangeMax = displayedRangeMax;
    }

    public Integer getDisplayedRangeIntervalCount() {
        return displayedRangeIntervalCount;
    }

    public void setDisplayedRangeIntervalCount(Integer displayedRangeIntervalCount) {
        this.displayedRangeIntervalCount = displayedRangeIntervalCount;
    }

    public String getDisplayedVariantType() {
        return displayedVariantType;
    }

    public void setDisplayedVariantType(String displayedVariantType) {
        this.displayedVariantType = displayedVariantType;
    }
}
