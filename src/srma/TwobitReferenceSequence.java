package srma;

import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;



public class TwobitReferenceSequence implements ReferenceSequenceFile {

	int index = 0;
	
	String chr;
	int size = 0;
	TwoBitGenomeReader tgr;
	SAMFileHeader samFileHeader;

	public TwobitReferenceSequence(TwoBitGenomeReader tgr, SAMFileHeader samFileHeader) {

		
		this.tgr = tgr;
		this.samFileHeader = samFileHeader;

	}

	public String getChr() {
		return chr;
	}

	@Override
	public void close() throws IOException {
		// do not do anyting
	}

	@Override
	public ReferenceSequence getSequence(String chr) {
		return null;// do not use
	}

	public ReferenceSequence getSequence(String chr, int start, int end) {

		try {
			int size = tgr.getReadSizes().get(chr);
			byte[] seq = new byte[size];
			if (start < 1) {
				start = 1;
			}
			if (end >= size - 1) {
				end = size - 1;
			}
			for (int n = start; n < end; n++) {

				seq[n-1] = (byte)tgr.getGenomeNuc(chr, n,true);

			}

			int index = samFileHeader.getSequenceIndex(chr);
			ReferenceSequence ref = new ReferenceSequence(chr, index, seq);
			return ref;

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	

	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		// TODO Auto-generated method stub
		return samFileHeader.getSequenceDictionary();
	}

	@Override
	public ReferenceSequence nextSequence() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean isIndexed() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start,
			long stop) {
		// TODO Auto-generated method stub
		return null;
	}

	

}
