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

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.PileupElementFilter;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

/**
 * GATK locus walker for computing a set of metrics for single nucleotide
 * variants (SNVs) from a given VCF file. These metrics are added to the
 * INFO column in the output VCF file and can be used for subsequent filtering
 * to obtain a higher precision set of variants.
 *
 * @author eldrid01
 */
@Requires(value = {DataSource.REFERENCE_ORDERED_DATA, DataSource.READS})
@Downsample(by = DownsampleType.NONE)
public class CalculateSNVMetrics extends LocusWalker<Integer, Integer>
{
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName = "countOverlappingReadPairs", shortName = "overlap", doc = "Handling of overlapping reads from the same fragment", required = false)
    private OverlappingReadPairPileupMode countOverlappingReadPairs = OverlappingReadPairPileupMode.DISCARD_IF_DISCORDANT_OR_USE_HIGHEST_BASE_QUALITY;

    @Argument(fullName="minimumBaseQuality", shortName="Q", doc="Minimum base quality at the variant position for reads to be included in the calculation of metrics.", required=false)
    private int minimumBaseQuality = 10;

    @Argument(fullName="minimumMappingQuality", shortName="M", doc="Minimum mapping quality for reads to be included in the calculation of metrics.", required=false)
    private int minimumMappingQuality = 1;

    @Argument(fullName="sample", doc="The sample(s) for which to compute SNV metrics. Can be specified multiple times for multiple samples; used to restrict the reads used in computing metrics, if not specified metrics are computed from all reads within the given BAM files.", required=false)
    private HashSet<String> samples;

    @Argument(fullName="controlSample", doc="The sample(s) for which to compute SNV metrics for the reference allele. Can be specified multiple times for multiple samples; if not specified metrics will be computed for all reference-supporting reads in the samples specified using the sample argument.", required=false)
    private HashSet<String> controlSamples;

    @Argument(fullName="excludeSupplementaryAlignments", doc="Whether to exclude supplementary alignments.", required=false)
    private boolean excludeSupplementaryAlignments;

    private static final int MAX_HOMOPOLYMER_LENGTH = 50;
    private static final int MAX_REPEAT_LENGTH = 30;

    private class BasePileupElementFilter implements PileupElementFilter
    {
    	private byte base;

    	public BasePileupElementFilter(byte base)
    	{
    		this.base = base;
    	}

    	@Override
		public boolean allow(PileupElement pileupElement)
    	{
    		return pileupElement.getBase() == base;
		}
    }

    /**
     * Filter for supplementary alignments.
     *
     * An alternative approach would be to use the --read_filter flag
     * although there doesn't appear to be a ReadFilter implementation
     * for excluding supplementary alignments.
     */
    private PileupElementFilter supplementaryAlignmentFilter = new PileupElementFilter() {
    	@Override
		public boolean allow(PileupElement pileupElement) {
    		return !excludeSupplementaryAlignments || !pileupElement.getRead().getSupplementaryAlignmentFlag();
		}
    };

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

    private ReferenceSequenceFile referenceSequence;

    /**
     * Initialization of count data structures, output header and read variants.
     */
    public void initialize()
    {
    	initializeSamples();
    	writeVCFHeader();
    	referenceSequence = getToolkit().getReferenceDataSource().getReference();
    }

    /**
     * Initialize samples.
     */
    private void initializeSamples()
    {
    	HashSet<String> allSamples = new HashSet<String>(getSampleDB().getSampleNames());

    	if (samples == null || samples.isEmpty())
    		samples = allSamples;
    	else
    		checkSamples(samples, allSamples);

    	if (controlSamples == null)
    		controlSamples = new HashSet<String>();
		checkSamples(controlSamples, allSamples);
    }

    /**
     * Check the given samples are among those in the set of all samples.
     *
     * @param samples
     * @param allSamples
     */
    private void checkSamples(HashSet<String> samples, HashSet<String> allSamples)
    {
        for (String sample : samples)
        {
        	if (!allSamples.contains(sample))
            {
        		String message = "Unrecognized sample: " + sample;
                logger.error(message);
                throw new RuntimeException(message);
            }
        }
    }

    /**
     * Initialize VCF header adding additional info lines for metrics.
     */
    private void writeVCFHeader()
    {
        List<String> rodName = Arrays.asList(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodName);

        Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);

        headerLines.add(new VCFInfoHeaderLine("Depth", 1, VCFHeaderLineType.Integer, "The number of reads covering the variant position including duplicates, supplementary records, reads that fall below minimum base and mapping quality thresholds, and from overlapping fragments."));
        headerLines.add(new VCFInfoHeaderLine("ReadCount", 1, VCFHeaderLineType.Integer, "The number of reads covering the variant position excluding duplicates, supplementary records and reads that fall below minimum base and mapping quality thresholds."));
        headerLines.add(new VCFInfoHeaderLine("VariantAlleleCount", 1, VCFHeaderLineType.Integer, "The variant allele count, i.e. the number of reads supporting the variant allele."));
        headerLines.add(new VCFInfoHeaderLine("VariantAlleleFrequency", 1, VCFHeaderLineType.Float, "The variant allele frequency, i.e. the fraction of reads supporting the variant allele."));

        if (!controlSamples.isEmpty())
        {
        	headerLines.add(new VCFInfoHeaderLine("ReadCountControl", 1, VCFHeaderLineType.Integer, "The number of reads covering the variant position in the control sample(s) excluding duplicates, supplementary records and reads that fall below minimum base and mapping quality thresholds."));
        	headerLines.add(new VCFInfoHeaderLine("VariantAlleleCountControl", 1, VCFHeaderLineType.Integer, "The variant allele count in the control sample(s)."));
        }

        headerLines.add(new VCFInfoHeaderLine("StrandBias", 1, VCFHeaderLineType.Float, "The strand bias for all reads covering the variant position."));
        headerLines.add(new VCFInfoHeaderLine("VariantStrandBias", 1, VCFHeaderLineType.Float, "The strand bias for variant-supporting reads."));
        headerLines.add(new VCFInfoHeaderLine("ReferenceStrandBias", 1, VCFHeaderLineType.Float, "The strand bias for reference-supporting reads."));

        headerLines.add(new VCFInfoHeaderLine("LowMapQual", 1, VCFHeaderLineType.Float, "The proportion of all reads from all samples at the variant position that have low mapping quality (less than " + minimumMappingQuality + ")."));

        headerLines.add(new VCFInfoHeaderLine("VariantBaseQual", 1, VCFHeaderLineType.Float, "The mean base quality at the variant position of variant reads."));
        headerLines.add(new VCFInfoHeaderLine("VariantBaseQualMedian", 1, VCFHeaderLineType.Float, "The median base quality at the variant position of variant reads."));

        headerLines.add(new VCFInfoHeaderLine("VariantMapQual", 1, VCFHeaderLineType.Float, "The mean mapping quality of variant reads."));
        headerLines.add(new VCFInfoHeaderLine("VariantMapQualMedian", 1, VCFHeaderLineType.Float, "The median mapping quality of variant reads."));
        headerLines.add(new VCFInfoHeaderLine("MapQualDiff", 1, VCFHeaderLineType.Float, "The difference in the mean mapping quality of variant and reference reads."));
        headerLines.add(new VCFInfoHeaderLine("MapQualDiffMedian", 1, VCFHeaderLineType.Float, "The difference in the median mapping quality of variant and reference reads."));

        headerLines.add(new VCFInfoHeaderLine("VariantMMQS", 1, VCFHeaderLineType.Float, "The mean mismatch quality sum for variant reads."));
        headerLines.add(new VCFInfoHeaderLine("VariantMMQSMedian", 1, VCFHeaderLineType.Float, "The median mismatch quality sum for variant reads."));
        headerLines.add(new VCFInfoHeaderLine("MMQSDiff", 1, VCFHeaderLineType.Float, "The difference in mean mismatch quality sum of variant and reference reads."));
        headerLines.add(new VCFInfoHeaderLine("MMQSDiffMedian", 1, VCFHeaderLineType.Float, "The difference in median mismatch quality sum of variant and reference reads."));

        headerLines.add(new VCFInfoHeaderLine("DistanceToAlignmentEnd", 1, VCFHeaderLineType.Float, "The mean shortest distance of the variant position within the read to either aligned end."));
        headerLines.add(new VCFInfoHeaderLine("DistanceToAlignmentEndMedian", 1, VCFHeaderLineType.Float, "The median shortest distance of the variant position within the read to either aligned end."));
        headerLines.add(new VCFInfoHeaderLine("DistanceToAlignmentEndMAD", 1, VCFHeaderLineType.Float, "The median absolute deviation of the shortest distance of the variant position within the read to either aligned end."));

        headerLines.add(new VCFInfoHeaderLine("HomopolymerLength", 1, VCFHeaderLineType.Integer, "The longest continuous homopolymer surrounding or adjacent to the variant position."));
        headerLines.add(new VCFInfoHeaderLine("Repeat", 1, VCFHeaderLineType.Integer, "The length of repetitive sequence adjacent to the variant position where repeats can be 1-, 2-, 3- or 4-mers."));

        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    private ReadBackedPileup filterSamples(ReadBackedPileup pileup, Collection<String> samplesToInclude)
    {
    	HashSet<String> samples = new HashSet<String>(pileup.getSamples());
    	samples.retainAll(samplesToInclude);
    	if (samples.isEmpty())
    		return new ReadBackedPileupImpl(pileup.getLocation());
    	else
    		return pileup.getPileupForSamples(samples);
    }

    @Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
	{
        String chromosome = context.getContig();
        long position = context.getPosition();
        byte referenceBase = ref.getBase();

        Collection<VariantContext> variants = tracker.getValues(variantCollection.variants, context.getLocation());
        if (variants.isEmpty()) return 0;

        for (VariantContext variant : variants)
        {
        	if (variant.isSNP())
        	{
        		ReadBackedPileup pileup = context.getBasePileup();
        		pileup = pileup.getFilteredPileup(supplementaryAlignmentFilter);

        		variant.getCommonInfo().putAttribute("Depth", filterSamples(pileup, samples).depthOfCoverage());

                int unfilteredReadCount = pileup.getBases().length;

                // apply mapping quality filter
                pileup = pileup.getMappingFilteredPileup(minimumMappingQuality);

                if (unfilteredReadCount > 0)
                {
                    variant.getCommonInfo().putAttribute("LowMapQual", 1.0 - (pileup.getBases().length / (double)unfilteredReadCount));
                }

                // apply base quality filter
                pileup = pileup.getBaseFilteredPileup(minimumBaseQuality);

                // apply filter for overlapping fragments
         		pileup = filterOverlappingFragements(pileup);

         		// apply filter for the specified sample(s)
                ReadBackedPileup samplePileup = filterSamples(pileup, samples);

                // collect reference metrics from control samples if these exist,
                // otherwise from reference-supporting reads in the same set of
                // samples used for computing the variant metrics
                ReadBackedPileup referencePileup = controlSamples.isEmpty() ? samplePileup : filterSamples(pileup, controlSamples);

                int readCount = samplePileup.getBases().length;
                variant.getCommonInfo().putAttribute("ReadCount", readCount);

                if (!variant.getReference().basesMatch(ref.getBases()))
                {
                	String message = "Mismatching reference base for variant " + getDisplayString(variant);
                	logger.error(message);
                    throw new RuntimeException(message);
                }
                int referenceBaseIndex = BaseUtils.simpleBaseToBaseIndex(referenceBase);

                if (variant.getAlternateAlleles().size() > 1)
                {
                    logger.warn("Warning multiple alternate alleles for variant " + getDisplayString(variant) + "; only computing metrics for allele with highest allele count.");
                }
                byte variantAllele = variant.getAltAlleleWithHighestAlleleCount().getBases()[0];
                int variantAlleleIndex = BaseUtils.simpleBaseToBaseIndex(variantAllele);

                int variantReadCount = samplePileup.getBaseCounts()[variantAlleleIndex];
                variant.getCommonInfo().putAttribute("VariantAlleleCount", variantReadCount);
                if (readCount > 0)
                {
                	variant.getCommonInfo().putAttribute("VariantAlleleFrequency", variantReadCount / (double)readCount);
                }

                int referenceReadCount = referencePileup.getBaseCounts()[referenceBaseIndex];

                if (!controlSamples.isEmpty())
                {
                    variant.getCommonInfo().putAttribute("ReadCountControl", referenceReadCount);
                    variant.getCommonInfo().putAttribute("VariantAlleleCountControl", referencePileup.getBaseCounts()[variantAlleleIndex]);
                }

                ReadBackedPileup negativeStrandPileup = samplePileup.getNegativeStrandPileup();

                if (readCount > 0)
                {
                    double strandBias = negativeStrandPileup.getBases().length / (double)readCount;
                    if (strandBias > 0.5) strandBias = 1.0 - strandBias;
                    variant.getCommonInfo().putAttribute("StrandBias", strandBias);
                }

                if (variantReadCount > 0)
                {
                    double variantStrandBias = negativeStrandPileup.getBaseCounts()[variantAlleleIndex] / (double)variantReadCount;
                    if (variantStrandBias > 0.5) variantStrandBias = 1.0 - variantStrandBias;
                    variant.getCommonInfo().putAttribute("VariantStrandBias", variantStrandBias);
                }

                if (referenceReadCount > 0)
                {
                    ReadBackedPileup referenceNegativeStrandPileup = referencePileup.getNegativeStrandPileup();
                    double referenceStrandBias = referenceNegativeStrandPileup.getBaseCounts()[referenceBaseIndex] / (double)referenceReadCount;
                    if (referenceStrandBias > 0.5) referenceStrandBias = 1.0 - referenceStrandBias;
                    variant.getCommonInfo().putAttribute("ReferenceStrandBias", referenceStrandBias);
                }

                DescriptiveStatistics variantBaseQualityStatistics = new DescriptiveStatistics();
                DescriptiveStatistics variantMappingQualityStatistics = new DescriptiveStatistics();
                DescriptiveStatistics variantMMQSStatistics = new DescriptiveStatistics();
                DescriptiveStatistics variantDistanceToAlignmentEndStatistics = new DescriptiveStatistics();

                ReadBackedPileup variantAllelePileup = samplePileup.getFilteredPileup(new BasePileupElementFilter(variantAllele));
                for (PileupElement pileupElement : variantAllelePileup)
                {
                	variantBaseQualityStatistics.addValue(pileupElement.getQual());
                	variantMappingQualityStatistics.addValue(pileupElement.getMappingQual());
                	variantMMQSStatistics.addValue(getMismatchQualitySum(pileupElement.getRead(), position));
                	variantDistanceToAlignmentEndStatistics.addValue(getDistanceToAlignmentEnd(chromosome, position, pileupElement));
                }

                if (variantBaseQualityStatistics.getN() > 0)
                {
                    variant.getCommonInfo().putAttribute("VariantBaseQual", variantBaseQualityStatistics.getMean());
                    variant.getCommonInfo().putAttribute("VariantBaseQualMedian", variantBaseQualityStatistics.getPercentile(50));
                }

                if (variantMappingQualityStatistics.getN() > 0)
                {
                    variant.getCommonInfo().putAttribute("VariantMapQual", variantMappingQualityStatistics.getMean());
                    variant.getCommonInfo().putAttribute("VariantMapQualMedian", variantMappingQualityStatistics.getPercentile(50));
                }

                DescriptiveStatistics referenceMappingQualityStatistics = new DescriptiveStatistics();
                DescriptiveStatistics referenceMMQSStatistics = new DescriptiveStatistics();

                ReadBackedPileup referenceBasePileup = referencePileup.getFilteredPileup(new BasePileupElementFilter(referenceBase));
                for (PileupElement pileupElement : referenceBasePileup)
                {
                	referenceMappingQualityStatistics.addValue(pileupElement.getMappingQual());
                	referenceMMQSStatistics.addValue(getMismatchQualitySum(pileupElement.getRead(), position));
                }

                if (variantMappingQualityStatistics.getN() > 0 && referenceMappingQualityStatistics.getN() > 0)
                {
                    variant.getCommonInfo().putAttribute("MapQualDiff", referenceMappingQualityStatistics.getMean() - variantMappingQualityStatistics.getMean());
                    variant.getCommonInfo().putAttribute("MapQualDiffMedian", referenceMappingQualityStatistics.getPercentile(50) - variantMappingQualityStatistics.getPercentile(50));
                }

                if (variantMMQSStatistics.getN() > 0)
                {
                    variant.getCommonInfo().putAttribute("VariantMMQS", variantMMQSStatistics.getMean());
                    variant.getCommonInfo().putAttribute("VariantMMQSMedian", variantMMQSStatistics.getPercentile(50));
                    if (referenceMMQSStatistics.getN() > 0)
                    {
                        variant.getCommonInfo().putAttribute("MMQSDiff", variantMMQSStatistics.getMean() - referenceMMQSStatistics.getMean());
                        variant.getCommonInfo().putAttribute("MMQSDiffMedian", variantMMQSStatistics.getPercentile(50) - referenceMMQSStatistics.getPercentile(50));
                    }
                }

                if (variantDistanceToAlignmentEndStatistics.getN() > 0)
                {
                    variant.getCommonInfo().putAttribute("DistanceToAlignmentEnd", variantDistanceToAlignmentEndStatistics.getMean());
                    variant.getCommonInfo().putAttribute("DistanceToAlignmentEndMedian", variantDistanceToAlignmentEndStatistics.getPercentile(50));
                    variant.getCommonInfo().putAttribute("DistanceToAlignmentEndMAD", getMedianAbsoluteDeviation(variantDistanceToAlignmentEndStatistics));
                }

                variant.getCommonInfo().putAttribute("HomopolymerLength", getHomopolymerLength(chromosome, position));
                variant.getCommonInfo().putAttribute("Repeat", getRepeatLength(chromosome, position));
        	}

        	vcfWriter.add(variant);
        }

        return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum)
	{
		return sum + value;
	}

	@Override
	public Integer reduceInit()
	{
		return 0;
	}

    public void onTraversalDone(Integer result)
    {
        logger.info("Number of locations traversed: " + result);
    }

    /**
     * Add a filter to the given read-backed pileup for overlapping read pairs.
     *
     * @param pileup the read-backed pileup.
     * @return
     */
    private ReadBackedPileup filterOverlappingFragements(ReadBackedPileup pileup)
    {
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

        return pileup;
    }

    /**
     * Returns a display string for the given variant.
     *
     * @param variantContext
     * @return
     */
    private String getDisplayString(VariantContext variant)
    {
        StringBuilder sb = new StringBuilder();
        sb.append(variant.getContig());
        sb.append(":");
        sb.append(variant.getStart());
        sb.append("_");
        sb.append(variant.getReference().getBaseString());
        for (Allele allele : variant.getAlternateAlleles())
        {
            sb.append("/");
            sb.append(allele.getBaseString());
        }
        return sb.toString();
    }

    /**
     * Calculates the median absolute deviation for the given statistics.
     *
     * @param statistics
     * @return
     */
    private double getMedianAbsoluteDeviation(DescriptiveStatistics statistics)
    {
        double median = statistics.getPercentile(50);
        DescriptiveStatistics absoluteDeviations = new DescriptiveStatistics();
        for (int i = 0; i < statistics.getN(); i++)
        {
            absoluteDeviations.addValue(Math.abs(statistics.getElement(i) - median));
        }
        return absoluteDeviations.getPercentile(50);
    }

    /**
     * Gets the mismatch quality sum for the given read alignment.
     *
     * @param record
     * @return
     */
    private int getMismatchQualitySum(SAMRecord record, long excludePosition)
    {
        byte[] readBases = record.getReadBases();
        byte[] baseQualities = record.getBaseQualities();

        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        byte[] referenceBases = referenceSequence.getSubsequenceAt(record.getReferenceName(), alignmentStart, alignmentEnd).getBases();

        int mmqs = 0;

        for (AlignmentBlock block : record.getAlignmentBlocks())
        {
            int readBlockStart = block.getReadStart() - 1;
            int referenceStart = block.getReferenceStart();
            int referenceBlockStart = referenceStart - alignmentStart;
            for (int i = 0; i < block.getLength(); i++)
            {
            	if ((referenceStart + i) == excludePosition) continue;
                byte referenceBase = referenceBases[referenceBlockStart + i];
                byte readBase = readBases[readBlockStart + i];
                if (!SequenceUtil.isNoCall(readBase) && !SequenceUtil.isNoCall(referenceBase) && !SequenceUtil.basesEqual(readBase, referenceBase))
                {
                    mmqs += baseQualities[readBlockStart + i];
                }
            }
        }

        return mmqs;
    }

    /**
     * Computes the distance between the variant position with a read and the
     * closest end of the aligned portion of the read.
     *
     * @param chromosome
     * @param position
     * @param pileupElement
     * @return
     */
    private int getDistanceToAlignmentEnd(String chromosome, long position, PileupElement pileupElement)
    {
    	GATKSAMRecord record = pileupElement.getRead();
    	int positionInRead = pileupElement.getOffset() + 1;
    	if (positionInRead != getPositionInRead(record, position))
    		logger.warn("Mismatching variant position within read");
        int alignedStartPositionInRead = getPositionInRead(record, record.getAlignmentStart());
        int alignedEndPositionInRead = getPositionInRead(record, record.getAlignmentEnd());
        int distance = Math.min(positionInRead - alignedStartPositionInRead, alignedEndPositionInRead - positionInRead) + 1;
        if (distance < 1)
        {
        	String message = "Error: unexpected distance to alignment start/end for read " + record.getReadName() + " at " + chromosome + ":" + position;
        	logger.error(message);
        	throw new RuntimeException(message);
        }
        return distance;
    }

    /**
     * Returns the position within the read that maps to the given position
     * in the reference sequence.
     *
     * @param record
     * @param referencePosition
     * @return
     */
    private int getPositionInRead(SAMRecord record, long referencePosition)
    {
        for (int i = 1; i <= record.getReadLength(); i++)
        {
            if (record.getReferencePositionAtReadPosition(i) == referencePosition) return i;
        }
        return 0;
    }

    /**
     * Returns the length of the longest homopolymer adjacent to or surrounding the given position.
     *
     * @param chromosome
     * @param position
     * @return
     */
    private int getHomopolymerLength(String chromosome, long position)
    {
        int referenceLength = referenceSequence.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        byte[] referenceBases = referenceSequence.getSubsequenceAt(chromosome, position, Math.min(referenceLength, position + MAX_HOMOPOLYMER_LENGTH)).getBases();
        byte referenceBase = referenceBases[0];
        byte rightBase = 0;
        int rightCount = 0;
        if (referenceBases.length > 1)
        {
            rightBase = referenceBases[1];
            rightCount = 1;
            for (int i = 2; i < referenceBases.length; i++)
            {
                if (!SequenceUtil.basesEqual(rightBase, referenceBases[i])) break;
                rightCount++;
            }
            if (SequenceUtil.basesEqual(rightBase, referenceBase)) rightCount++;
        }

        referenceBases = referenceSequence.getSubsequenceAt(chromosome, Math.max(1, position - MAX_HOMOPOLYMER_LENGTH), position).getBases();
        byte leftBase = 0;
        int leftCount = 0;
        if (referenceBases.length > 1)
        {
            leftBase = referenceBases[referenceBases.length - 2];
            leftCount = 1;
            for (int i = referenceBases.length - 3; i >= 0; i--)
            {
                if (!SequenceUtil.basesEqual(leftBase, referenceBases[i])) break;
                leftCount++;
            }
            if (SequenceUtil.basesEqual(leftBase, referenceBase)) leftCount++;
        }

        if (leftCount > 0 && rightCount > 0 && SequenceUtil.basesEqual(leftBase, referenceBase) && SequenceUtil.basesEqual(rightBase, referenceBase))
        {
            return leftCount + rightCount - 1;
        }
        else
        {
            return Math.max(leftCount, rightCount);
        }
    }

    /**
     * Returns the length of repetitive sequence adjacent to the variant
     * position where repeats can be 1-, 2-, 3- or 4-mers.
     *
     * @param chromosome
     * @param position
     * @return
     */
    private int getRepeatLength(String chromosome, long position)
    {
        int repeat1 = getRepeatCount(1, chromosome, position);
        int repeat2 = getRepeatCount(2, chromosome, position);
        int repeat3 = getRepeatCount(3, chromosome, position);
        int repeat4 = getRepeatCount(4, chromosome, position);

        int repeat = 0;
        if (repeat1 > 1 || repeat2 > 1 || repeat3 > 1 || repeat4 > 1)
        {
            repeat = Math.max(repeat1,  repeat2 * 2);
            repeat = Math.max(repeat, repeat3 * 3);
            repeat = Math.max(repeat, repeat4 * 4);
        }

        return repeat;
    }

    /**
     * Returns the maximum number of repeats of the specified length adjacent to
     * the given position.
     *
     * @param n the repeat length, e.g. n = 2 for dinucleotide repeats
     * @param chromosome
     * @param position
     * @return
     */
    private int getRepeatCount(int n, String chromosome, long position)
    {
        int referenceLength = referenceSequence.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        int count = 0;

        if (position < referenceLength)
        {
            byte[] referenceBases = referenceSequence.getSubsequenceAt(chromosome, position + 1, Math.min(referenceLength, position + MAX_REPEAT_LENGTH)).getBases();
            count = Math.max(count, getRepeatCount(n, referenceBases));
        }

        if (position > 1)
        {
            byte[] referenceBases = referenceSequence.getSubsequenceAt(chromosome, Math.max(1, position - MAX_REPEAT_LENGTH), position - 1).getBases();
            ArrayUtils.reverse(referenceBases);
            count = Math.max(count, getRepeatCount(n, referenceBases));
        }

        return count;
    }

    /**
     * Returns the number of repeats of length n in the given sequence.
     *
     * @param n the repeat length, e.g. n = 2 for dinucleotide repeats
     * @param bases
     * @return
     */
    private int getRepeatCount(int n, byte[] bases)
    {
        if (bases.length < n) return 0;
        int count = 1;
        while (true)
        {
            int pos = count * n;
            if ((pos + n) > bases.length) return count;
            for (int i = 0; i < n; i++, pos++)
            {
                if (bases[i] != bases[pos]) return count;
            }
            count++;
        }
    }

}
