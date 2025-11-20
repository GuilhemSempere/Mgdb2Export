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
package htsjdk.variant.variantcontext.writer;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * This class is pretty much a clone of HTSJDK's
 * org.broadinstitute.sting.utils.variantcontext.writer.VCFWriter. It was added
 * in order to be able to force the default locale (en_US) to be used when
 * formatting numbers (methods formatVCFDouble & formatQualValue)
 *
 * @author SEMPERE
 */
public class CustomVCFWriter extends IndexingVariantContextWriter {

    /**
     * Instantiates a new custom vcf writer.
     *
     * @param location the location
     * @param output the output
     * @param refDict the ref dict
     * @param enableOnTheFlyIndexing the enable on the fly indexing
     * @param doNotWriteGenotypes the do not write genotypes
     * @param allowMissingFieldsInHeader the allow missing fields in header
     */
    public CustomVCFWriter(File location, OutputStream output, SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes, boolean allowMissingFieldsInHeader) {
        super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        mHeader = null;
        intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();
        this.doNotWriteGenotypes = doNotWriteGenotypes;
        this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
    }

    /* (non-Javadoc)
     * @see htsjdk.variant.variantcontext.writer.IndexingVariantContextWriter#writeHeader(htsjdk.variant.vcf.VCFHeader)
     */
    public void writeHeader(VCFHeader header) {
        mHeader = writeHeader(header, ((Writer) (mWriter)), doNotWriteGenotypes, getVersionLine(), getStreamName());
    }
    
    public void accountForHeaderWithoutWritingIt(VCFHeader header) {
        mHeader = writeHeader(header, new StringWriter(), doNotWriteGenotypes, getVersionLine(), getStreamName());
    }

    /**
     * Gets the version line.
     *
     * @return the version line
     */
    public static final String getVersionLine() {
        return VERSION_LINE;
    }

    /**
     * Write header.
     *
     * @param header the header
     * @param writer the writer
     * @param doNotWriteGenotypes the do not write genotypes
     * @param versionLine the version line
     * @param streamNameForError the stream name for error
     * @return the VCF header
     */
    public static VCFHeader writeHeader(VCFHeader header, Writer writer, boolean doNotWriteGenotypes, String versionLine, String streamNameForError) {
        header = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInSortedOrder()) : header;
        try {
            writer.write((new StringBuilder()).append(versionLine).append("\n").toString());
            Iterator i$ = header.getMetaDataInSortedOrder().iterator();
            do {
                if (!i$.hasNext()) {
                    break;
                }
                VCFHeaderLine line = (VCFHeaderLine) i$.next();
                if (!VCFHeaderVersion.isFormatString(line.getKey())) {
                    writer.write("##");
                    writer.write(line.toString());
                    writer.write("\n");
                }
            } while (true);
            writer.write("#");
            boolean isFirst = true;
            htsjdk.variant.vcf.VCFHeader.HEADER_FIELDS field;
            for (i$ = header.getHeaderFields().iterator(); i$.hasNext(); writer.write(field.toString())) {
                field = (htsjdk.variant.vcf.VCFHeader.HEADER_FIELDS) i$.next();
                if (isFirst) {
                    isFirst = false;
                } else {
                    writer.write("\t");
                }
            }

            if (header.hasGenotypingData()) {
                writer.write("\t");
                writer.write("FORMAT");
                String sample;
                for (i$ = header.getGenotypeSamples().iterator(); i$.hasNext(); writer.write(sample)) {
                    sample = (String) i$.next();
                    writer.write("\t");
                }

            }
            writer.write("\n");
            writer.flush();
        } catch (IOException e) {
            throw new Error((new StringBuilder()).append("IOException writing the VCF header to ").append(streamNameForError).toString(), e);
        }
        return header;
    }

    /* (non-Javadoc)
     * @see htsjdk.variant.variantcontext.writer.IndexingVariantContextWriter#close()
     */
    public void close() {
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new Error((new StringBuilder()).append("Unable to close ").append(getStreamName()).toString(), e);
        }
        super.close();
    }

    /* (non-Javadoc)
     * @see htsjdk.variant.variantcontext.writer.IndexingVariantContextWriter#add(htsjdk.variant.variantcontext.VariantContext)
     */
    public void add(VariantContext vc) {
        if (mHeader == null) {
            throw new IllegalStateException((new StringBuilder()).append("The VCF Header must be written before records can be added: ").append(getStreamName()).toString());
        }
        if (doNotWriteGenotypes) {
            vc = (new VariantContextBuilder(vc)).noGenotypes().make();
        }
        try {
            super.add(vc);
            Map alleleMap = buildAlleleMap(vc);
            mWriter.write(vc.getChr());
            mWriter.write("\t");
            mWriter.write(String.valueOf(vc.getStart()));
            mWriter.write("\t");
            String ID = vc.getID();
            mWriter.write(ID);
            mWriter.write("\t");
            String refString = vc.getReference().getDisplayString();
            mWriter.write(refString);
            mWriter.write("\t");
            if (vc.isVariant()) {
                Allele altAllele = vc.getAlternateAllele(0);
                String alt = altAllele.getDisplayString();
                mWriter.write(alt);
                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alt = altAllele.getDisplayString();
                    mWriter.write(",");
                    mWriter.write(alt);
                }

            } else {
                mWriter.write(".");
            }
            mWriter.write("\t");
            if (!vc.hasLog10PError()) {
                mWriter.write(".");
            } else {
                mWriter.write(formatQualValue(vc.getPhredScaledQual()));
            }
            mWriter.write("\t");
            String filters = getFilterString(vc);
            mWriter.write(filters);
            mWriter.write("\t");
            Map infoFields = new TreeMap();
            Iterator i$ = vc.getAttributes().entrySet().iterator();
            do {
                if (!i$.hasNext()) {
                    break;
                }
                java.util.Map.Entry field = (java.util.Map.Entry) i$.next();
                String key = (String) field.getKey();
                if (!mHeader.hasInfoLine(key)) {
                    fieldIsMissingFromHeaderError(vc, key, "INFO");
                }
                String outputValue = formatVCFField(field.getValue());
                if (outputValue != null) {
                    infoFields.put(key, outputValue);
                }
            } while (true);
            writeInfoString(infoFields);
            GenotypesContext gc = vc.getGenotypes();
            if (gc.isLazyWithData() && (((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String)) {
                mWriter.write("\t");
                mWriter.write(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
            } else {
                List genotypeAttributeKeys = calcVCFGenotypeKeys(vc, mHeader);
                if (!genotypeAttributeKeys.isEmpty()) {
                    i$ = genotypeAttributeKeys.iterator();
                    do {
                        if (!i$.hasNext()) {
                            break;
                        }
                        String format = (String) i$.next();
                        if (!mHeader.hasFormatLine(format)) {
                            fieldIsMissingFromHeaderError(vc, format, "FORMAT");
                        }
                    } while (true);
                    String genotypeFormatString = ParsingUtils.join(":", genotypeAttributeKeys);
                    mWriter.write("\t");
                    mWriter.write(genotypeFormatString);
                    addGenotypeData(vc, alleleMap, genotypeAttributeKeys);
                }
            }
            mWriter.write("\n");
            mWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException((new StringBuilder()).append("Unable to write the VCF object to ").append(getStreamName()).toString(), e);
        }
    }

    /**
     * Builds the allele map.
     *
     * @param vc the vc
     * @return the map
     */
    private static Map buildAlleleMap(VariantContext vc) {
        Map alleleMap = new HashMap(vc.getAlleles().size() + 1);
        alleleMap.put(Allele.NO_CALL, ".");
        List alleles = vc.getAlleles();
        for (int i = 0; i < alleles.size(); i++) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

    /**
     * Gets the filter string.
     *
     * @param vc the vc
     * @return the filter string
     */
    private final String getFilterString(VariantContext vc) {
        if (vc.isFiltered()) {
            Iterator i$ = vc.getFilters().iterator();
            do {
                if (!i$.hasNext()) {
                    break;
                }
                String filter = (String) i$.next();
                if (!mHeader.hasFilterLine(filter)) {
                    fieldIsMissingFromHeaderError(vc, filter, "FILTER");
                }
            } while (true);
            return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
        }
        if (vc.filtersWereApplied()) {
            return "PASS";
        } else {
            return ".";
        }
    }

    /**
     * Format qual value.
     *
     * @param qual the qual
     * @return the string
     */
    private String formatQualValue(double qual) {
        String s = String.format(VCFConstants.VCF_LOCALE, "%.2f", new Object[]{
            Double.valueOf(qual)
        });
        if (s.endsWith(".00")) {
            s = s.substring(0, s.length() - ".00".length());
        }
        return s;
    }

    /**
     * Write info string.
     *
     * @param infoFields the info fields
     * @throws IOException Signals that an I/O exception has occurred.
     */
    private void writeInfoString(Map infoFields)
            throws IOException {
        if (infoFields.isEmpty()) {
            mWriter.write(".");
            return;
        }
        boolean isFirst = true;
        Iterator i$ = infoFields.entrySet().iterator();
        do {
            if (!i$.hasNext()) {
                break;
            }
            java.util.Map.Entry entry = (java.util.Map.Entry) i$.next();
            if (isFirst) {
                isFirst = false;
            } else {
                mWriter.write(";");
            }
            String key = (String) entry.getKey();
            mWriter.write(key);
            if (!((String) entry.getValue()).equals("")) {
                VCFInfoHeaderLine metaData = mHeader.getInfoHeaderLine(key);
                if (metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0) {
                    mWriter.write("=");
                    mWriter.write(/*String.format(VCFConstants.VCF_LOCALE, */(String) entry.getValue()/*)*/);
                }
            }
        } while (true);
    }

    /**
     * Adds the genotype data.
     *
     * @param vc the vc
     * @param alleleMap the allele map
     * @param genotypeFormatKeys the genotype format keys
     * @throws IOException Signals that an I/O exception has occurred.
     */
    private void addGenotypeData(VariantContext vc, Map alleleMap, List genotypeFormatKeys)
            throws IOException {
        int ploidy = vc.getMaxPloidy(2);
        for (Iterator j$ = mHeader.getGenotypeSamples().iterator(); j$.hasNext();) {
            String sample = (String) j$.next();
            mWriter.write("\t");
            Genotype g = vc.getGenotype(sample);
            if (g == null) {
                g = GenotypeBuilder.createMissing(sample, ploidy);
            }
            List attrs = new ArrayList(genotypeFormatKeys.size());
            Iterator i$ = genotypeFormatKeys.iterator();
            do {
                if (!i$.hasNext()) {
                    break;
                }
                String field = (String) i$.next();
                if (field.equals("GT")) {
                    if (!g.isAvailable()) {
                        throw new Error("GTs cannot be missing for some samples if they are available for others in the record");
                    }
                    writeAllele(g.getAllele(0), alleleMap);
                    int i = 1;
                    while (i < g.getPloidy()) {
                        mWriter.write(g.isPhased() ? "|" : "/");
                        writeAllele(g.getAllele(i), alleleMap);
                        i++;
                    }
                } else {
                    String outputValue;
                    if (field.equals("FT")) {
                        outputValue = g.isFiltered() ? g.getFilters() : "PASS";
                    } else {
                        IntGenotypeFieldAccessors.Accessor accessor = intGenotypeFieldAccessors.getAccessor(field);
                        if (accessor != null) {
                            int intValues[] = accessor.getValues(g);
                            if (intValues == null) {
                                outputValue = ".";
                            } else if (intValues.length == 1) {
                                outputValue = Integer.toString(intValues[0]);
                            } else {
                                StringBuilder sb = new StringBuilder();
                                sb.append(intValues[0]);
                                for (int i = 1; i < intValues.length; i++) {
                                    sb.append(",");
                                    sb.append(intValues[i]);
                                }

                                outputValue = sb.toString();
                            }
                        } else {
                            Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : ".";
                            VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(field);
                            if (metaData != null) {
                                int numInFormatField = Math.min(10000, metaData.getCount(vc));	// limit applied to avoid OutOfMemoryError on some MIXED type variants where the value ended up being huge (this change does not seem to have any effect on the output)
                                if (numInFormatField > 1 && ".".equals(val == null ? "." : val)) {
                                    StringBuilder sb = new StringBuilder(".");
                                    for (int i = 1; i < numInFormatField; i++)
                                        sb.append(",.");

                                    val = sb.toString();
                                }
                            }
                            outputValue = formatVCFField(val);
                        }
                    }
                    if (outputValue != null) {
                        attrs.add(outputValue);
                    }
                }
            } while (true);
            
            int i;
            
            // removing this optimization which suppresses the last attribute value if it is a missing value (we prefer to be more strictly compliant with VCF format)
//            for (i = attrs.size() - 1; i >= 0 && isMissingValue((String) attrs.get(i)); i--) {
//                attrs.remove(i);
//            }

            i = 0;
            while (i < attrs.size()) {
                if (i > 0 || genotypeFormatKeys.contains("GT")) {
                    mWriter.write(":");
                }
                mWriter.write((String) attrs.get(i));
                i++;
            }
        }
    }

    /**
     * Checks if is missing value.
     *
     * @param s the s
     * @return true, if is missing value
     */
    private boolean isMissingValue(String s) {
        return countOccurrences(".".charAt(0), s) + countOccurrences(',', s) == s.length();
    }

    /**
     * Write allele.
     *
     * @param allele the allele
     * @param alleleMap the allele map
     * @throws IOException Signals that an I/O exception has occurred.
     */
    private void writeAllele(Allele allele, Map alleleMap)
            throws IOException {
        String encoding = (String) alleleMap.get(allele);
        if (encoding == null) {
            throw new TribbleException.InternalCodecException((new StringBuilder()).append("Allele ").append(allele).append(" is not an allele in the variant context").toString());
        } else {
            mWriter.write(encoding);
            return;
        }
    }

    /**
     * Format vcf double.
     *
     * @param d the d
     * @return the string
     */
    public static final String formatVCFDouble(double d) {
        String format;
        if (d < 1.0D) {
            if (d < 0.01D) {
                if (Math.abs(d) >= 9.9999999999999995E-21D) {
                    format = "%.3e";
                } else {
                    return "0.00";
                }
            } else {
                format = "%.3f";
            }
        } else {
            format = "%.2f";
        }
        return String.format(VCFConstants.VCF_LOCALE, format, new Object[]{Double.valueOf(d)});
    }

    /**
     * Format vcf field.
     *
     * @param val the val
     * @return the string
     */
    public static String formatVCFField(Object val) {
        String result;
        if (val == null) {
            result = ".";
        } else if (val instanceof Double) {
            result = formatVCFDouble(((Double) val).doubleValue());
        } else if (val instanceof Boolean) {
            result = ((Boolean) val).booleanValue() ? "" : null;
        } else if (val instanceof List) {
            result = formatVCFField(((Object) (((List) val).toArray())));
        } else if (val.getClass().isArray()) {
            int length = Array.getLength(val);
            if (length == 0) {
                return formatVCFField(null);
            }
            StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
            for (int i = 1; i < length; i++) {
                sb.append(",");
                sb.append(formatVCFField(Array.get(val, i)));
            }

            result = sb.toString();
        } else {
            result = val.toString();
        }
        return result;
    }

    /**
     * Calc vcf genotype keys.
     *
     * @param vc the vc
     * @param header the header
     * @return the list
     */
    public static List calcVCFGenotypeKeys(VariantContext vc, VCFHeader header) {
        Set keys = new HashSet();
        boolean sawGoodGT = false;
        boolean sawGoodQual = false;
        boolean sawGenotypeFilter = false;
        boolean sawDP = false;
        boolean sawAD = false;
        boolean sawPL = false;
        Iterator i$ = vc.getGenotypes().iterator();
        do {
            if (!i$.hasNext()) {
                break;
            }
            Genotype g = (Genotype) i$.next();
            keys.addAll(g.getExtendedAttributes().keySet());
            if (g.isAvailable()) {
                sawGoodGT = true;
            }
            if (g.hasGQ()) {
                sawGoodQual = true;
            }
            if (g.hasDP()) {
                sawDP = true;
            }
            if (g.hasAD()) {
                sawAD = true;
            }
            if (g.hasPL()) {
                sawPL = true;
            }
            if (g.isFiltered()) {
                sawGenotypeFilter = true;
            }
        } while (true);
        if (sawGoodQual) {
            keys.add("GQ");
        }
        if (sawDP) {
            keys.add("DP");
        }
        if (sawAD) {
            keys.add("AD");
        }
        if (sawPL) {
            keys.add("PL");
        }
        if (sawGenotypeFilter) {
            keys.add("FT");
        }
        List sortedList = ParsingUtils.sortList(new ArrayList(keys));
        if (sawGoodGT) {
            List newList = new ArrayList(sortedList.size() + 1);
            newList.add("GT");
            newList.addAll(sortedList);
            sortedList = newList;
        }
        if (sortedList.isEmpty() && header.hasGenotypingData()) {
            return Collections.singletonList("GT");
        } else {
            return sortedList;
        }
    }

    /**
     * Count occurrences.
     *
     * @param c the char
     * @param s the string
     * @return the occurence count
     */
    private static int countOccurrences(char c, String s) {
        int count = 0;
        for (int i = 0; i < s.length(); i++) {
            count += s.charAt(i) != c ? 0 : 1;
        }

        return count;
    }

    /**
     * Field is missing from header error.
     *
     * @param vc the vc
     * @param id the id
     * @param field the field
     */
    private void fieldIsMissingFromHeaderError(VariantContext vc, String id, String field) {
        if (!allowMissingFieldsInHeader) {
            throw new RuntimeException((new StringBuilder()).append("Malformed VCF Header: Key ").append(id).append(" found in VariantContext field ").append(field).append(" at ").append(vc.getChr()).append(":").append(vc.getStart()).append(" but this key isn't defined in the VCFHeader. The GATK now requires all VCFs to have").append(" complete VCF headers by default. This error can be disabled with the engine argument").append(" -U LENIENT_VCF_PROCESSING").toString());
        } else {

        }
    }

    /**
     * The Constant VERSION_LINE.
     */
    private static final String VERSION_LINE;

    /**
     * The writer.
     */
    protected final BufferedWriter mWriter = new BufferedWriter(new OutputStreamWriter(getOutputStream()));

    /**
     * The do not write genotypes.
     */
    protected final boolean doNotWriteGenotypes;

    /**
     * The vcf header.
     */
    protected VCFHeader mHeader;

    /**
     * The allow missing fields in header.
     */
    private final boolean allowMissingFieldsInHeader;

    /**
     * The int genotype field accessors.
     */
    private IntGenotypeFieldAccessors intGenotypeFieldAccessors;

    /**
     * The Constant QUAL_FORMAT_STRING.
     */
    private static final String QUAL_FORMAT_STRING = "%.2f";

    /**
     * The Constant QUAL_FORMAT_EXTENSION_TO_TRIM.
     */
    private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

    static {
        VERSION_LINE = (new StringBuilder()).append("##").append(VCFHeaderVersion.VCF4_1.getFormatString()).append("=").append(VCFHeaderVersion.VCF4_1.getVersionString()).toString();
    }

	@Override
	public void setHeader(VCFHeader header) {
		this.mHeader = header;
	}
}
