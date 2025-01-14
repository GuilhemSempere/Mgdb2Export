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
public class MgdbVcfFieldPlotRequest extends MgdbDensityRequest {

    private String vcfField;
    
    public MgdbVcfFieldPlotRequest(){
        super();
    }

	public String getVcfField() {
		return vcfField;
	}

	public void setVcfField(String vcfField) {
		this.vcfField = vcfField;
	}
}