package org.cruk.gatk.walkers;

import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.util.SequenceUtil;

/**
 * GATK read walker for adding UMI tags to records within a BAM file where the
 * unique molecular index or barcode is the first N bases in the read sequence.
 *
 * This utility was written to add 6-base barcodes from the beginning of each
 * read to the alignment records so that these can be taken into consideration
 * when duplicate marking with Picard MarkDuplicates.
 *
 * Note that this does not add UMI tags to secondary or supplementary alignment
 * records.
 *
 * @author eldrid01
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class AddUMITags extends ReadWalker<GATKSAMRecord, SAMFileWriter> implements NanoSchedulable
{
    @Output(doc="Write output to this BAM filename instead of STDOUT")
    GATKSAMFileWriter out;

    @Argument(fullName="umiLength", doc="The number of bases at the beginning of the read that consitute the UMI or barcode.", required=true)
    private int umiLength;

    @Argument(fullName="umiTag", doc="The tag to use for the UMI or barcode in the output BAM file.", required=false)
    private String umiTag = "BX";

    public void initialize()
    {
        GenomeAnalysisEngine toolkit = getToolkit();
        if (toolkit != null)
        {
            SAMFileHeader outputHeader = toolkit.getSAMFileHeader().clone();
            out.writeHeader(outputHeader);
            out.setPresorted(true);
        }
    }

    @Override
    public GATKSAMRecord map(ReferenceContext referenceContext, GATKSAMRecord record, RefMetaDataTracker metaDataTracker)
    {
        if (!record.isSecondaryOrSupplementary())
        {
            String sequence = record.getReadString();
            if (record.getReadNegativeStrandFlag()) sequence = SequenceUtil.reverseComplement(sequence);
            String umi = sequence.substring(0, umiLength);
            record.setAttribute(umiTag, umi);
        }
        return record;
    }

    @Override
    public SAMFileWriter reduce(GATKSAMRecord record, SAMFileWriter output)
    {
        output.addAlignment(record);
        return output;
    }

    @Override
    public SAMFileWriter reduceInit()
    {
        return out;
    }
}
