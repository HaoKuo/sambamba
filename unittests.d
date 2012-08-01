import samheader;
import bamfile;
import samfile;
import sam.recordparser;
import bgzfrange;

import validation.samheader;
import validation.alignment;

import utils.samheadermerger;

import sam.serialize;

import std.path;
import std.stdio;
import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.math;

unittest {

    writeln("Testing extracting SAM header...");
    auto fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    auto bf = BamFile(fn);
    assert(bf.header.format_version == "1.3");
    assert(bf.header.sorting_order == SortingOrder.coordinate);
    assert(bf.header.sequences.length == 2);
    assert(bf.header.getSequenceIndex("chr1") == 0);
    assert(bf.header.sequences["chr2"].length == 1584);

    fn = buildPath(dirName(__FILE__), "test", "data", "bins.bam");
    bf = BamFile(fn);
    assert(bf.header.sorting_order == SortingOrder.unknown);
    assert(bf.header.sequences.length == 3);
    assert(bf.header.read_groups.length == 0);
    assert(bf.header.getSequenceIndex("large") == 2);
    assert(bf.header.sequences["small"].length == 65536);

    writeln("Testing alignment parsing...");
    fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    bf = BamFile(fn);
    auto alignments = bf.alignments;

    {
    auto read = alignments.front;
    assert(equal(read.sequence, "CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA"));
    assert(equal(map!"cast(char)(a + 33)"(read.phred_base_quality),
                "<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;"));
    assert(bf.reference(read.ref_id).name == "chr1");
    assert(read.read_name == "EAS56_57:6:190:289:82");
    assert(read.flag == 69);
    assert(read.position == 99);
    assert(read.mapping_quality == 0);
    alignments.popFront();
    alignments.popFront();
    assert(alignments.front.cigarString == "35M");
    assert(toSam(alignments.front, bf.reference_sequences) == "EAS51_64:3:190:727:308	99	chr1	103	99	35M	=	263	195	GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG	<<<<<<<<<<<<<<<<<<<<<<<<<<<::<<<844	MF:i:18	Aq:i:73	NM:i:0	UQ:i:0	H0:i:1	H1:i:0");
    assert(bf.header.getSequenceIndex("chr1") == read.ref_id);
    }

    assert(bf.alignments.front.read_name == "EAS56_57:6:190:289:82");

    writeln("Testing tag parsing...");
    fn = buildPath(dirName(__FILE__), "test", "data", "tags.bam");
    bf = BamFile(fn);
    foreach (alignment; bf.alignments) {
        auto read_name = alignment.read_name;
        assert(read_name[0..4] == "tag_");
        char[] tag;
        read_name = read_name[4..$];
        while (read_name[0] != ':') {
            tag ~= read_name[0];
            read_name = read_name[1..$];
        }
        read_name = read_name[1..$];
        string value = toSam(alignment.tags[tag.idup]);
        if (read_name != value) {
            writeln("tag: ", tag, "\tread_name: ", read_name, "\tvalue: ", value);
            writeln("value bam_typeid: ", alignment.tags[tag.idup].bam_typeid);
		}

        assert(read_name == value);
    }

    writeln("Testing exception handling...");
    fn = buildPath(dirName(__FILE__), "test", "data", "duplicated_block_size.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "no_block_size.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "wrong_extra_gzip_length.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "wrong_bc_subfield_length.bam");
    assertThrown!BgzfException(reduce!"a+b.sequence_length"(0,BamFile(fn).alignments));
    fn = buildPath(dirName(__FILE__), "test", "data", "corrupted_zlib_archive.bam");
    assertThrown!ZlibException(BamFile(fn));

    writeln("Testing random access...");
    fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    bf = BamFile(fn);

    void compareWithNaiveApproach(int beg, int end) {

        auto refseq = array(bf["chr1"][beg .. end]);

        auto naive = array(filter!((Alignment a) { 
                         return bf.reference(a.ref_id).name == "chr1" &&
                                a.position < end &&
                                a.position + a.basesCovered() > beg; })
                            (bf.alignments));
        if (!equal(naive, refseq)) {
            writeln(beg);
            writeln(end);
            writeln(array(map!"a.read_name"(refseq)));
            writeln(array(map!"a.read_name"(naive)));
        }
        assert(equal(refseq, naive));
    }

    compareWithNaiveApproach(1400, 1500);
    compareWithNaiveApproach(  10,  123);
    compareWithNaiveApproach( 135, 1236);
    compareWithNaiveApproach(1350, 3612);
    compareWithNaiveApproach( 643, 1732);
    compareWithNaiveApproach( 267, 1463);
    compareWithNaiveApproach(   0,   30);
    compareWithNaiveApproach(1363, 1612);
    compareWithNaiveApproach( 361, 1231);
    compareWithNaiveApproach( 322,  612);
    compareWithNaiveApproach( 912,  938);
    compareWithNaiveApproach(-100, 3000);
    compareWithNaiveApproach(-100,  100);
    compareWithNaiveApproach(-100, 1000);
    compareWithNaiveApproach(   0, 1900);
    compareWithNaiveApproach(   1,  279);

    writeln("Testing Value code...");
    Value v = 5;
    assert(v.is_integer);
    assert(toSam(v) == "i:5");
    assert(v == 5);
    assert(v == "5");
    assert(v != [1,2,3]);
    v = "abc";
    assert(v.is_string);
    assert(toSam(v) == "Z:abc");
    assert(v == "abc");
    v = [1, 2, 3];
    assert(v.is_numeric_array);
    assert(toSam(v) == "B:i,1,2,3");
    assert(v == [1,2,3]);
    assert(v == "[1, 2, 3]");
    v = [1.5, 2.3, 17.0];
    assert(v.is_numeric_array);
    assert(toSam(v) == "B:f,1.5,2.3,17");
    assert(approxEqual(to!(float[])(v), [1.5, 2.3, 17]));
    v = 5.6;
    assert(v.is_float);
    assert(toSam(v) == "f:5.6");
    assert(approxEqual(to!float(v), 5.6));
    v = -17;
    assert(v.is_signed);
    assert(toSam(v) == "i:-17");
    assert(v == -17);
    assert(v == "-17");
    v = 297u;
    assert(v.is_unsigned);
    assert(toSam(v) == "i:297");
    assert(v == 297);
    assert(v == "297");

    short[] array_of_shorts = [4, 5, 6];
    v = array_of_shorts;
    assert(v.is_numeric_array);
    assert(toSam(v) == "B:s,4,5,6");
    assert(to!(short[])(v) == array_of_shorts);
    assert(v == [4,5,6]);
    assert(v == "[4, 5, 6]");

    v = null;
    assert(v.is_nothing);

    v = "0eabcf123";
    v.setHexadecimalFlag();
    assert(v.is_hexadecimal_string);    
    assert(v == "0eabcf123");

    writeln("Test parseAlignmentLine/toSam functions...");
    fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    bf = BamFile(fn);
    foreach (read; bf.alignments) {
        auto line = toSam(read, bf.reference_sequences);
        auto read2 = parseAlignmentLine(line, bf.header);
        if (read != read2) {
            writeln(read.read_name);
        }
        assert(read == read2);
    }

    fn = buildPath(dirName(__FILE__), "test", "data", "tags.bam");
    bf = BamFile(fn);
    foreach (read; bf.alignments) {
        auto line = toSam(read, bf.reference_sequences);
        auto read2 = parseAlignmentLine(line, bf.header);
        if (read != read2 && isValid(bf.header, read)) {
            writeln(read.read_name);
        }
        assert(read == read2 || !isValid(bf.header, read));
    }

}

void main() {
}
