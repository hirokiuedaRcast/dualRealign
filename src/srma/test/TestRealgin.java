/*
 * Copyright Hiroki Ueda
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


package srma.test;

import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.realignimpl.Realignment;

public class TestRealgin {
	
	

	
	public static void main(String[] arg){
		

		

		String normalbamf = "/home/spark/todaitoptest/testdata0308/LUAD-311-N_TDv3_genome.bam";
		String tumorbamf = "/home/spark/todaitoptest/testdata0308/LUAD-311-T_TDv3_genome.bam";
		
		
		String twobitref = "/home/spark/todaitoptest/ref/hg38.2bit";
		
		String targetRegion = "/home/spark/todaitoptest/ref/S3035822_Covered_sm_lift_hg38.bed";
		
		
		String outdir = "/home/spark/todaitoptest/testdata0308/realign";
		
		List<String> l = new ArrayList<String>();
		
		add(l,"-n",normalbamf);
		add(l,"-t",tumorbamf);
		add(l,"-r",twobitref);
		
		add(l,"-ct",targetRegion);
	
		add(l,"-o",outdir);
		
		add(l,"-nt","4");
		
				
		
		String[] ar = l.toArray(new String[l.size()]);
		
		Realignment.main(ar);
		
	}

	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

}
