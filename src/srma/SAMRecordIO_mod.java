package srma;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;


public class SAMRecordIO_mod 
{
   
	List<SAMRecord> inputs;
	List<SAMRecord> outputs;
	int idx = 0;
    public SAMRecordIO_mod(List<SAMRecord> inputs,SAMSequenceDictionary referenceDictionary)
        throws Exception
    {
     
    	this.inputs = inputs;
    	this.outputs = new ArrayList<SAMRecord>();
    	idx= 0;
    }

	public void output(AlignRecord rec) {
		
		rec.record.setAttribute("YU","toSRMA");
		outputs.add(rec.record);
		
	}

	public boolean hasNextAlignRecord() {
		
		return idx >= 0 && idx<inputs.size();
	}

	public AlignRecord getNextAlignRecord() {
		
		SAMRecord sam = inputs.get(idx);
		idx++;
		return new AlignRecord(sam,null,0);
		
	}
    
}
