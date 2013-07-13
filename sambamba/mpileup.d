/*
    This file is part of Sambamba.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.mpileup;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.snpcallers.maq;

import bio.core.format;

import std.stdio;

version(standalone) {
    int main(string[] args) {
        return mpileup_main(args);
    }
}

void writeInfoField(SNP)(scope void delegate(const(char)[]) sink, auto ref SNP snp) {
    sink.writeGenotype(snp);
    sink.writeGenotypeLikelihoods(snp);
    sink.writeGenotypeLikelihoodsOfHeterogeneousPloidy(snp);
    sink.writeRoundedGenotypeLikelihoods(snp);
    sink.writeGenotypePosteriorProbabilities(snp);
    sink.writeConditionalGenotypeQuality(snp);
    sink.writeHaplotypeQualities(snp);
    sink.writePhaseSet(snp);
    sink.writePhasingQuality(snp);
    sink.writeExpectedAlternateAlleleCounts(snp);
    sink.writeRootMeanSquareMappingQuality(snp);
}

void toVCF(SNP)(scope void delegate(const(char)[]) sink, auto ref SNP snp) {
    // CHROM
    sink.write(snp.sample);
    sink.write('\t');

    // POS
    sink.write(snp.position);

    // ID
    sink.write("\t.\t");

    // REF
    sink.write(snp.reference_base);
    sink.write('\t');

    // ALT
    if (snp.genotype.base1 != snp.reference_base)
        sink.write(snp.genotype.base1);

    if (snp.genotype.is_heterozygous &&
        snp.genotype.base2 != snp.reference_base)
    {
        sink.write(snp.genotype.base2);
    }

    // QUAL
    sink.write('\t');
    sink.write(snp.quality);

    // FILTER
    sink.write("\t.\t");

    // INFO
    // TODO sink.writeInfoField(snp);

    sink.write('\n');
}

void printVCFHeader(Writer)(ref Writer writer) {
    // TODO
}

void printVCFLine(Writer, SNP)(ref Writer writer, auto ref SNP snp) {
    toVCF((const(char)[] s) { writer.put(s); }, snp);
}

int mpileup_main(string[] args) {
    try {
        auto bam = new BamReader(args[1]);
        auto writer = stdout.lockingTextWriter;
        
        foreach (pileup; pileupChunks(bam.reads, useMD, 128_000_000)) {
            auto caller = new MaqSnpCaller();
            caller.minimum_call_quality = 20.0f;
            caller.minimum_base_quality = 15;
            foreach (snp; caller.findSNPs(pileup)) {
                writer.printVCFLine(snp);
            }
        }

    } catch (Throwable e) {
        stderr.writeln(e.msg);
    }
}
