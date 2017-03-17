package jp.ac.utokyo.rcast.realignimpl;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Indel {

	int pos;
	public int len;
	public boolean insersion;
	
	int depth =1;
	
	int normal =0;
	int tumor = 0;
	int total =0;
	
	public Indel(SAMRecord sam) {
		
		int n= sam.getAlignmentStart();
		for (CigarElement ce : sam.getCigar().getCigarElements()) {
			
			
			
			if (ce.getOperator().equals(CigarOperator.D)) {
				pos = n;
				len = ce.getLength();
				break;
			}
			if (ce.getOperator().equals(CigarOperator.I)) {
				pos = n;
				len = ce.getLength();
				break;
			}
			
			if(ce.getOperator().consumesReferenceBases()){
				pos = pos + ce.getLength();
			}
				
 
		}
		
	}

	public void incTotal() {
		
		depth++;
		
	}

	public void inc(boolean normalb) {
		
		if(normalb){
			normal++;
		}else{
			tumor++;
		}
		total++;
	}

	public int getCount() {
		
		return total;
	}

	public double getAlelleFrequency() {
		
		return (double)total/(double)depth;
		
	}

	public boolean observedBoth() {
		
		return (normal >2) && (tumor >2);
	}

}
