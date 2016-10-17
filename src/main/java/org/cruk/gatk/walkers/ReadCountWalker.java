/*
 * MIT License
 * 
 * Copyright (c) 2016 CRUK Cambridge Institute - Bioinformatics Core
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.cruk.gatk.walkers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.PileupElementFilter;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

/**
 * GATK locus walker for getting read counts for each base, optionally writing
 * reference and variant-supporting read counts for a set of SNVs given in a VCF file.
 *
 * @author eldrid01
 */
@Downsample(by = DownsampleType.NONE)
public class ReadCountWalker extends LocusWalker<Integer, Long> implements PileupElementFilter
{
    @Output(fullName="output_file", shortName="o", doc="Write read counts to this file instead of STDOUT")
    PrintStream out;

    @Argument(fullName="sample", doc="The sample(s) to output read counts for (can be specified multiple times for multiple samples)", required=false)
    private List<String> samples;

    @Argument(fullName="variantsFile", doc="File containing SNVs in VCF format; if specified will output alternate allele(s) and optionally the reference and alternate allele counts and allele fraction (see outputRefAltCounts and outputAlleleFraction options).", required=false)
    private File variantsFile;

    @Argument(fullName = "countOverlappingReadPairs", shortName = "overlap", doc = "Handling of overlapping reads from the same fragment", required = false)
    private OverlappingReadPairPileupMode countOverlappingReadPairs = OverlappingReadPairPileupMode.RETAIN_BOTH;

    @Argument(fullName="properPairsOnly", doc="Only include reads mapped in proper pair.", required=false)
    private boolean properPairsOnly = false;

    @Argument(fullName="minimumBaseQuality", shortName="Q", doc="Minimum base quality to be included in read counts.", required=false)
    private int minimumBaseQuality = 0;

    @Argument(fullName="minimumMappingQuality", shortName="M", doc="Minimum mapping quality for reads to be included in counts.", required=false)
    private int minimumMappingQuality = 0;

    @Argument(fullName="outputDepth", doc="Output the depth, i.e. total count for all bases.", required=false)
    private boolean outputDepth = false;

    @Argument(fullName="outputUnfilteredDepth", doc="Output the unfiltered depth, i.e. total count for all bases prior to filtering based on base and mapping qualities, proper read pairs, etc.", required=false)
    private boolean outputUnfilteredDepth = false;

    @Argument(fullName="outputSeparateStrandCounts", doc="Output separate read counts for forward and reverse strands.", required=false)
    private boolean outputSeparateStrandCounts = false;

    @Argument(fullName="outputACGTCounts", doc="Output base counts for A, C, G and T bases; only applies when variant file supplied, otherwise these are output automatically.", required=false)
    private boolean outputACGTCounts = false;

    @Argument(fullName="outputNCounts", doc="Output counts for N bases.", required=false)
    private boolean outputNCounts = false;

    @Argument(fullName="outputRefAltCounts", doc="Output base counts for reference base and alternate allele(s); only applies when variant file supplied.", required=false)
    private boolean outputRefAltCounts = false;

    @Argument(fullName="outputAlleleFraction", doc="Output allele fraction; only applies when variant file supplied.", required=false)
    private boolean outputAlleleFraction = false;

    @Argument(fullName="minimumDepthPerSample", doc="Minimum depth per sample specified as comma-separated key=value pairs, e.g. SampleA=6,SampleB=10,...; note that all samples need to meet the minimum depth for the read counts to be output.", required=false)
    private String minimumDepthPerSample;

    @Argument(fullName="defaultMinimumDepthPerSample", doc="Default value for the minimum depth for each sample if not explicitly specified using minimumDepthPerSample.", required=false)
    private int defaultMinimumDepthPerSample = 0;

    @Argument(fullName="maximumDepthPerSample", doc="Maximum depth per sample specified as comma-separated key=value pairs, e.g. SampleA=100,SampleB=150,...; note that if any sample exceeds the maximum depth for that sample the read counts will not be output for that position.", required=false)
    private String maximumDepthPerSample;

    @Argument(fullName="defaultMaximumDepthPerSample", doc="Default value for the maximum depth for each sample if not explicitly specified using maximumDepthPerSample.", required=false)
    private int defaultMaximumDepthPerSample = Integer.MAX_VALUE;

    @Argument(fullName="ignoreReferenceSequenceWhenAssessingVariants", doc="Ignore the reference sequence when assessing variants using minimum variant count thresholds, e.g. when looking for heterozygosity within a single sample; if set to true, the secondmost frequentbase is taken as the variant allele.")
    private boolean ignoreReferenceSequenceWhenAssessingVariants = false;

    @Argument(fullName="minimumVariantCounts", doc="Minimum number of reads required to support a variant allele in each sample specified as comma-separated key/value pairs, e.g. SampleA=3,SampleB=4,...; note that only one of the samples needs to meet the variant allele count threshold for the read counts to be output for all samples.", required=false)
    private String minimumVariantCounts;

    @Argument(fullName="defaultMinimumVariantCount", doc="Default value for the minimum number of reads required to support variant allele for any sample where this is not explicitly specified using minimumVariantCounts.", required=false)
    private int defaultMinimumVariantCount = 0;

    @Argument(fullName="minimumVariantCountAcrossSamples", doc="Minimum number of reads required to support a variant allele when comparing across all samples; the counts for each base are summed for all specified samples and the read counts for each sample are output if the minor (variant) allele count is at or above this value.", required=false)
    private int minimumVariantCountAcrossSamples = 0;

    private static final String[] COUNT_COLUMNS = new String[] { "A", "C", "G", "T", "N", "Ref", "Alt" };
    private static final Map<Character, Integer> baseIndex = new HashMap<Character, Integer>();
    static
    {
        baseIndex.put('A', 0);
        baseIndex.put('C', 1);
        baseIndex.put('G', 2);
        baseIndex.put('T', 3);
        baseIndex.put('N', 4);
    }

    private Map<String, Integer> depthLookup = new HashMap<String, Integer>();
    private Map<String, Integer> unfilteredDepthLookup = new HashMap<String, Integer>();

    private Map<String, int[]> forwardCountLookup = new HashMap<String, int[]>();
    private Map<String, int[]> reverseCountLookup = new HashMap<String, int[]>();

    private boolean outputVariants = false;
    private Map<String, Map<Long, Variant[]>> variantLookup = new HashMap<String, Map<Long, Variant[]>>();

    private Map<String, Integer> minimumDepthLookup;
    private Map<String, Integer> maximumDepthLookup;
    private Map<String, Integer> minimumVariantCountLookup;
    private boolean depthFiltering = false;
    private boolean variantCountFiltering = false;

    /**
     * Enumeration of options for handling overlapping read pairs
     * at the locus of interest.
     */
    public enum OverlappingReadPairPileupMode
    {
        /**
         * Count both reads of a pair from the same fragment.
         */
        RETAIN_BOTH,
        /**
         * Discard both reads if they have different base calls at the locus of interest,
         * or use read with highest mapping quality if they have the same base.
         */
        DISCARD_IF_DISCORDANT_OR_USE_HIGHEST_MAPPING_QUALITY,
        /**
         * Discard both reads if they have different base calls at the locus of interest,
         * or use read with highest base quality if they have the same base.
         */
        DISCARD_IF_DISCORDANT_OR_USE_HIGHEST_BASE_QUALITY,
        /**
         * Use read with highest mapping quality.
         */
        RETAIN_HIGHEST_MAPPING_QUALITY,
        /**
         * Use read with highest base quality.
         */
        RETAIN_HIGHEST_BASE_QUALITY
    }

    /**
     * Initialization of count data structures, output header and read variants.
     */
    public void initialize()
    {
        initializeOutputOptions();
        initializeSampleCountLookups();
        initializeSampleDepthVariantCountLookups();
        printHeader();
        readVariants();
    }

    /**
     * Initialize default values for various output options if these haven't been
     * explicitly set.
     */
    private void initializeOutputOptions()
    {
        outputVariants = variantsFile != null;
        if (outputVariants)
        {
            if (!outputRefAltCounts && !outputAlleleFraction)
            {
                outputRefAltCounts = true;
                outputAlleleFraction = true;
            }
        }
        else
        {
            outputACGTCounts = true;
            outputRefAltCounts = false;
            outputAlleleFraction = false;
        }
    }

    /**
     * Initialize sample count lookups.
     */
    private void initializeSampleCountLookups()
    {
        Set<String> sampleSet = new HashSet<String>();
        List<String> sampleList = new ArrayList<String>();

        if (samples != null)
        {
            for (String sample : samples)
            {
                if (sample == null || sample.isEmpty()) continue;
                if (getSampleDB().getSample(sample) == null)
                {
                    String message = "Unrecognized sample: " + sample;
                	logger.error(message);
                    throw new RuntimeException(message);
                }
                if (sampleSet.contains(sample)) continue;
                sampleSet.add(sample);
                sampleList.add(sample);
            }
        }

        if (sampleList.isEmpty())
        {
            sampleList.addAll(getSampleDB().getSampleNames());
        }

        samples = sampleList;

        for (String sample : samples)
        {
            int n = baseIndex.size();
            forwardCountLookup.put(sample, new int[n]);
            reverseCountLookup.put(sample, new int[n]);
        }
    }

    /**
     * Print header.
     */
    private void printHeader()
    {
        out.print("Chromosome\tPosition\tReference base");
        if (outputVariants) out.print("\tAlternate non-reference allele");

        for (String sample : samples)
        {
            if (outputUnfilteredDepth)
            {
                out.print("\tDepth unfiltered");
                if (samples.size() > 1) out.print(" " + sample);
            }

            if (outputDepth)
            {
                out.print("\tDepth");
                if (samples.size() > 1) out.print(" " + sample);
            }

            int n = baseIndex.size();
            boolean[] include = new boolean[COUNT_COLUMNS.length];
            for (int i = 0; i < include.length; i++) include[i] = false;
            if (outputACGTCounts)
            {
                for (int i = 0; i < n; i++) include[i] = true;
            }
            include[n - 1] = outputNCounts;
            if (outputRefAltCounts)
            {
                include[n] = true;
                include[n + 1] = true;
            }

            for (int i = 0; i < COUNT_COLUMNS.length; i++)
            {
                if (include[i])
                {
                    if (outputSeparateStrandCounts)
                    {
                        out.print("\t" + COUNT_COLUMNS[i] + " forward count");
                        if (samples.size() > 1) out.print(" " + sample);
                        out.print("\t" + COUNT_COLUMNS[i] + " reverse count");
                        if (samples.size() > 1) out.print(" " + sample);
                    }
                    else
                    {
                        out.print("\t" + COUNT_COLUMNS[i] + " count");
                        if (samples.size() > 1) out.print(" " + sample);
                    }
                }
            }

            if (outputAlleleFraction)
            {
                out.print("\tAllele fraction");
                if (samples.size() > 1) out.print(" " + sample);
            }
        }

        out.println();
    }

    /**
     * Initialize sample-based depth and variant count lookups.
     */
    private void initializeSampleDepthVariantCountLookups()
    {
        minimumDepthLookup = parseSampleValues("minimumDepthPerSample", minimumDepthPerSample, defaultMinimumDepthPerSample);
        maximumDepthLookup = parseSampleValues("maximumDepthPerSample", maximumDepthPerSample, defaultMaximumDepthPerSample);
        minimumVariantCountLookup = parseSampleValues("minimumVariantCounts", minimumVariantCounts, defaultMinimumVariantCount);
        depthFiltering = minimumDepthPerSample != null || maximumDepthPerSample != null || defaultMinimumDepthPerSample != 0 || defaultMaximumDepthPerSample != Integer.MAX_VALUE;
        variantCountFiltering = minimumVariantCounts != null || defaultMinimumVariantCount > 0;
    }

    /**
     * Parses arguments given for sample-based integer values of the form
     * SM1=5,SM2=3,SM3=10,... and creates a lookup.
     *
     * @param name the argument name.
     * @param values
     * @param defaultValue
     * @return
     */
    private Map<String, Integer> parseSampleValues(String name, String values, int defaultValue)
    {
        Map<String, Integer> lookup = new HashMap<String, Integer>();

        for (String sample : samples)
        {
            lookup.put(sample, defaultValue);
        }

        if (values == null || values.isEmpty()) return lookup;

        String[] tokens = values.split("[,;]");
        for (String token : tokens)
        {
            String[] subtokens = token.split("=");
            if (subtokens.length != 2)
            {
                String message = "Incorrect format for " + name + " parameter: " + values;
            	logger.error(message);
                throw new RuntimeException(message);
            }
            String sample = subtokens[0];
            if (!lookup.containsKey(sample))
            {
                String message = "Unrecognized sample " + sample + " in setting for " + name + ": " + values;
            	logger.error(message);
                throw new RuntimeException(message);
            }
            try
            {
                Integer value = Integer.parseInt(subtokens[1]);
                lookup.put(sample, value);
            }
            catch (NumberFormatException e)
            {
                String message = "Incorrect format for " + name + " parameter (expecting integer values): " + values;
            	logger.error(message);
                throw new RuntimeException(message);
            }
        }

        return lookup;
    }

    /**
     * Checks if the current position passes the sample-based depth filter.
     *
     * Note that only one of the samples needs to be in the acceptable range.
     *
     * @return <code>true</code> if the depth of any one sample is within the acceptable range.
     */
    private boolean passDepthFilter()
    {
        if (depthFiltering)
        {
            for (String sample : samples)
            {
                int depth = depthLookup.get(sample);
                if (depth < minimumDepthLookup.get(sample) || depth > maximumDepthLookup.get(sample)) return false;
            }
        }
        return true;
    }

    /**
     * Checks if the current position passes sample-based minimum variant count
     * filter.
     *
     * Note that only one of the samples needs to be heterozygous with a minor
     * allele count greater than or equal to the specified threshold.
     *
     * @return <code>true</code> if the variant allele count of any one sample is at or above the specified threshold.
     */
    private boolean passVariantCountFilter(char referenceBase)
    {
        if (variantCountFiltering)
        {
            int n = baseIndex.size() - 1;
            int referenceBaseIndex = baseIndex.get(referenceBase);

            for (String sample : samples)
            {
                int[] forwardCounts = forwardCountLookup.get(sample);
                int[] reverseCounts = reverseCountLookup.get(sample);

                int minimumVariantCount = minimumVariantCountLookup.get(sample);
                int aboveMinimumVariantCount = 0;

                for (int i = 0; i < n; i++)
                {
                    if (!ignoreReferenceSequenceWhenAssessingVariants && i == referenceBaseIndex) continue;
                    int count = forwardCounts[i] + reverseCounts[i];
                    if (count >= minimumVariantCount) aboveMinimumVariantCount++;
                }

                if (ignoreReferenceSequenceWhenAssessingVariants)
                {
                    if (aboveMinimumVariantCount > 1) return true;
                }
                else
                {
                    if (aboveMinimumVariantCount > 0) return true;
                }
            }

            return false;
        }
        else
        {
            return true;
        }
    }

    /**
     * Checks if the current position passes cross-sample minimum variant count
     * filter.
     *
     * Combines read counts for all samples and then determines if the minor
     * allele count is at or above the specified threshold.
     *
     * @param referenceBase the reference base at the current position
     * @return <code>true</code> if the variant allele count across all samples is at or above the specified threshold.
     */
    private boolean passVariantCountAcrossSamplesFilter(char referenceBase)
    {
        if (minimumVariantCountAcrossSamples == 0) return true;

        int n = baseIndex.size() - 1;
        int[] counts = new int[n];

        for (String sample : samples)
        {
            int[] forwardCounts = forwardCountLookup.get(sample);
            int[] reverseCounts = reverseCountLookup.get(sample);

            for (int i = 0; i < n; i++)
            {
                counts[i] += forwardCounts[i] + reverseCounts[i];
            }
        }

        int referenceBaseIndex = baseIndex.get(referenceBase);
        int aboveMinimumVariantCount = 0;

        for (int i = 0; i < n; i++)
        {
            if (!ignoreReferenceSequenceWhenAssessingVariants && i == referenceBaseIndex) continue;
            if (counts[i] >= minimumVariantCountAcrossSamples) aboveMinimumVariantCount++;
        }

        if (ignoreReferenceSequenceWhenAssessingVariants)
        {
            return aboveMinimumVariantCount > 1;
        }
        else
        {
            return aboveMinimumVariantCount > 0;
        }
    }

    @Override
    public boolean allow(PileupElement pileupElement)
    {
        // filtering for reads in proper pair if specified
        return !properPairsOnly || pileupElement.getRead().getProperPairFlag();
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
        String chromosome = context.getContig();
        long position = context.getPosition();
        char referenceBase = (char)ref.getBase();

        ReadBackedPileup pileup = context.getBasePileup();
        for (String sample : samples)
        {
            unfilteredDepthLookup.put(sample, 0);
            ReadBackedPileup samplePileup = pileup.getPileupForSample(sample);
            if (samplePileup != null)
            {
                unfilteredDepthLookup.put(sample, samplePileup.getBases().length);
            }
        }

        // handling of overlapping read pairs
        if (countOverlappingReadPairs == OverlappingReadPairPileupMode.DISCARD_IF_DISCORDANT_OR_USE_HIGHEST_BASE_QUALITY)
        {
            pileup = pileup.getOverlappingFragmentFilteredPileup(true, true);
        }
        else if (countOverlappingReadPairs == OverlappingReadPairPileupMode.DISCARD_IF_DISCORDANT_OR_USE_HIGHEST_MAPPING_QUALITY)
        {
            pileup = pileup.getOverlappingFragmentFilteredPileup(true, false);
        }
        else if (countOverlappingReadPairs == OverlappingReadPairPileupMode.RETAIN_HIGHEST_BASE_QUALITY)
        {
            pileup = pileup.getOverlappingFragmentFilteredPileup(false, true);
        }
        else if (countOverlappingReadPairs == OverlappingReadPairPileupMode.RETAIN_HIGHEST_MAPPING_QUALITY)
        {
            pileup = pileup.getOverlappingFragmentFilteredPileup(false, false);
        }

        pileup = pileup.getFilteredPileup(this).getMappingFilteredPileup(minimumMappingQuality).getBaseFilteredPileup(minimumBaseQuality);

        for (String sample : samples)
        {
            depthLookup.put(sample, 0);

            int[] forwardCounts = forwardCountLookup.get(sample);
            int[] reverseCounts = reverseCountLookup.get(sample);
            for (int i = 0; i < baseIndex.size(); i++)
            {
                forwardCounts[i] = 0;
                reverseCounts[i] = 0;
            }

            ReadBackedPileup samplePileup = pileup.getPileupForSample(sample);
            if (samplePileup != null)
            {
                depthLookup.put(sample, samplePileup.getBases().length);

                Iterator<PileupElement> pileupIterator = samplePileup.iterator();
                while (pileupIterator.hasNext())
                {
                    PileupElement pileupElement = pileupIterator.next();

                    char base = (char)pileupElement.getBase();
                    int index = baseIndex.get(base);

                    boolean negativeStrand = pileupElement.getRead().getReadNegativeStrandFlag();
                    if (negativeStrand)
                    {
                        reverseCounts[index]++;
                    }
                    else
                    {
                        forwardCounts[index]++;
                    }
                }
            }
        }

        if (!passDepthFilter()) return 1;

        if (variantCountFiltering && minimumVariantCountAcrossSamples > 0)
        {
            // both variant count tests in play - only require one to pass
            if (!(passVariantCountFilter(referenceBase) || passVariantCountAcrossSamplesFilter(referenceBase))) return 1;
        }
        else
        {
            // only one or neither of the variant count tests in play - require both to return true
            if (!passVariantCountFilter(referenceBase)) return 1;
            if (!passVariantCountAcrossSamplesFilter(referenceBase)) return 1;
        }

        if (outputVariants)
        {
            Variant[] variants = getVariants(chromosome, position);
            if (variants == null)
            {
                outputVariantReadCounts(chromosome, position, referenceBase, null);
            }
            else
            {
                for (Variant variant : variants)
                {
                    if (referenceBase != variant.getReferenceBase())
                    {
                        String message = "Reference base at " + chromosome + ":" + position + " differs in the reference sequence (" + referenceBase + ") and variant file (" + variant.getReferenceBase() + ")";
                    	logger.error(message);
                        throw new RuntimeException(message);
                    }
                    outputVariantReadCounts(chromosome, position, referenceBase, variant);
                }
            }
        }
        else
        {
            outputReadCounts(chromosome, position, referenceBase);
        }

        return 1;
    }

    @Override
    public Long reduceInit()
    {
        return 0l;
    }

    @Override
    public Long reduce(Integer value, Long sum)
    {
        return value + sum;
    }

    public void onTraversalDone(Long c)
    {
        logger.info("Number of locations traversed: " + c);
    }

    /**
     * Output read counts for the given position.
     *
     * @param chromosome
     * @param position
     * @param referenceBase
     */
    private void outputReadCounts(String chromosome, long position, char referenceBase)
    {
        out.print(chromosome + "\t" + position + "\t" + referenceBase);
        for (String sample : samples)
        {
            if (outputUnfilteredDepth) out.print("\t" + unfilteredDepthLookup.get(sample));
            if (outputDepth) out.print("\t" + depthLookup.get(sample));

            if (outputACGTCounts || outputNCounts)
            {
                int[] forwardCounts = forwardCountLookup.get(sample);
                int[] reverseCounts = reverseCountLookup.get(sample);

                int n = baseIndex.size();
                int start = outputACGTCounts ? 0 : (n - 1);
                int end = outputNCounts ? n : (n - 1);

                for (int i = start; i < end; i++)
                {
                    if (outputSeparateStrandCounts)
                    {
                        out.print("\t" + forwardCounts[i]);
                        out.print("\t" + reverseCounts[i]);
                    }
                    else
                    {
                        out.print("\t" + (forwardCounts[i] + reverseCounts[i]));
                    }
                }
            }
        }

        out.println();
    }

    /**
     * Output read counts for the given variant.
     *
     * @param chromosome
     * @param position
     * @param referenceBase
     * @param variant
     */
    private void outputVariantReadCounts(String chromosome, long position, char referenceBase, Variant variant)
    {
        out.print(chromosome + "\t" + position + "\t" + referenceBase);

        out.print("\t");
        if (variant != null)
        {
            boolean first = true;
            for (char alternateAllele : variant.getAlternateAlleles())
            {
                if (first)
                {
                    first = false;
                }
                else
                {
                    out.print(",");
                }
                out.print(alternateAllele);
            }
        }

        for (String sample : samples)
        {
            if (outputUnfilteredDepth) out.print("\t" + unfilteredDepthLookup.get(sample));

            int depth = depthLookup.get(sample);
            if (outputDepth) out.print("\t" + depthLookup.get(sample));

            int[] forwardCounts = forwardCountLookup.get(sample);
            int[] reverseCounts = reverseCountLookup.get(sample);

            if (outputACGTCounts || outputNCounts)
            {
                int n = baseIndex.size();
                int start = outputACGTCounts ? 0 : (n - 1);
                int end = outputNCounts ? n : (n - 1);

                for (int i = start; i < end; i++)
                {
                    if (outputSeparateStrandCounts)
                    {
                        out.print("\t" + forwardCounts[i]);
                        out.print("\t" + reverseCounts[i]);
                    }
                    else
                    {
                        out.print("\t" + (forwardCounts[i] + reverseCounts[i]));
                    }
                }
            }

            if (variant == null)
            {
                if (outputRefAltCounts)
                {
                    out.print("\t\t");
                    if (outputSeparateStrandCounts) out.print("\t\t");
                }
                if (outputAlleleFraction)
                {
                    out.print("\t");
                }
            }
            else
            {
                int refForwardCount = forwardCounts[baseIndex.get(referenceBase)];
                int refReverseCount = reverseCounts[baseIndex.get(referenceBase)];
                int refCount = refForwardCount + refReverseCount;
                int altForwardCount = 0;
                int altReverseCount = 0;
                for (char alternateAllele : variant.getAlternateAlleles())
                {
                    if (baseIndex.containsKey(alternateAllele))
                    {
                        int index = baseIndex.get(alternateAllele);
                        altForwardCount += forwardCounts[index];
                        altReverseCount += reverseCounts[index];
                    }
                }
                int altCount = altForwardCount + altReverseCount;

                if (outputRefAltCounts)
                {
                    if (outputSeparateStrandCounts)
                    {
                        out.print("\t" + refForwardCount + "\t" + refReverseCount);
                        out.print("\t" + altForwardCount + "\t" + altReverseCount);
                    }
                    else
                    {
                        out.print("\t" + refCount);
                        out.print("\t" + altCount);
                    }
                }

                if (outputAlleleFraction)
                {
                    double alleleFraction = depth == 0 ? 0.0 : (double)altCount / depth;
                    out.printf("\t%.3f", alleleFraction);
                }
            }
        }

        out.println();
    }

    /**
     * Reads variants from the VCF file.
     */
    private void readVariants()
    {
        if (variantsFile == null) return;
        BufferedReader reader = null;
        try
        {
            reader = new BufferedReader(new FileReader(variantsFile));
            String line = null;
            while ((line = reader.readLine()) != null)
            {
                if (line.startsWith("#")) continue;

                String[] fields = line.split("\t");

                String chromosome = fields[0];
                long position = Long.parseLong(fields[1]);

                if (fields[3].length() != 1)
                {
                    throw new IllegalArgumentException("Unexpected reference base at line: " + line);
                }
                char referenceBase = fields[3].charAt(0);

                String[] tokens = fields[4].split(",");
                int n = tokens.length;
                char[] alternateAlleles = new char[n];
                for (int i = 0; i < n; i++)
                {
                    if (tokens[i].length() != 1)
                    {
                        throw new IllegalArgumentException("Unexpected variant base at line: " + line);
                    }
                    alternateAlleles[i] = tokens[i].charAt(0);
                }

                Variant variant = new Variant(referenceBase, alternateAlleles);
                addVariantToLookup(chromosome, position, variant);
            }
        }
        catch (IOException e)
        {
            logger.fatal("Error reading variants file");
            throw new RuntimeException(e);
        }
        finally
        {
            try
            {
                if (reader != null) reader.close();
            }
            catch (IOException e)
            {
            }
        }
    }

    /**
     * Adds the given variant to the position-based lookup.
     *
     * @param chromosome
     * @param position
     * @param referenceBase
     * @param alternateAlleles
     */
    private void addVariantToLookup(String chromosome, long position, Variant variant)
    {
        Map<Long, Variant[]> lookup = variantLookup.get(chromosome);
        if (lookup == null)
        {
            lookup = new HashMap<Long, Variant[]>();
            variantLookup.put(chromosome, lookup);
        }

        Variant[] variants = lookup.get(position);
        if (variants == null)
        {
            variants = new Variant[] { variant };
        }
        else
        {
            // check if the same variant has already been added to the lookup
            for (Variant v : variants)
            {
                if (variant.equals(v)) return;
            }
            Variant[] newVariants = new Variant[variants.length + 1];
            System.arraycopy(variants, 0, newVariants, 0, variants.length);
            newVariants[variants.length] = variant;
            variants = newVariants;
        }
        lookup.put(position, variants);
    }

    /**
     * Returns the set of variants at the given chromosome/position.
     *
     * @param chromosome
     * @param position
     * @return
     */
    private Variant[] getVariants(String chromosome, long position)
    {
        if (variantsFile == null) return null;

        Map<Long, Variant[]> lookup = variantLookup.get(chromosome);
        if (lookup == null) return null;

        Variant[] variants = lookup.get(position);
        if (variants == null) return null;

        return variants;
    }

    /**
     * Class representing a variant with reference base and alternate alleles.
     *
     * Chromosome and position information contained in the lookup and not duplicated
     * here to save on space.
     *
     * @author eldrid01
     */
    private class Variant
    {
        private char referenceBase;
        private char[] alternateAlleles;

        public Variant(char referenceBase, char[] alternateAlleles)
        {
            this.referenceBase = referenceBase;
            this.alternateAlleles = alternateAlleles;
        }

        public char getReferenceBase()
        {
            return referenceBase;
        }

        public char[] getAlternateAlleles()
        {
            return alternateAlleles;
        }

        @Override
        public boolean equals(Object obj)
        {
            if (obj instanceof Variant)
            {
                Variant other = (Variant)obj;
                if (referenceBase != other.referenceBase) return false;
                if (alternateAlleles.length != other.alternateAlleles.length) return false;
                for (int i = 0; i < alternateAlleles.length; i++)
                {
                    if (alternateAlleles[i] != other.alternateAlleles[i]) return false;
                }
                return true;
            }
            else
            {
                return false;
            }
        }

    }
}
