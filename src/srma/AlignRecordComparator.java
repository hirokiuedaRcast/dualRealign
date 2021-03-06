package srma;

import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.Comparator;


public class AlignRecordComparator implements Comparator<AlignRecord> 
{
    private SAMRecordComparator comp = null;

    public AlignRecordComparator()
    {
        comp = new SAMRecordCoordinateComparator();
    }

    public int compare(AlignRecord o1, AlignRecord o2)
    {
        return this.comp.compare(o1.record, o2.record);
    }

    public boolean equals(AlignRecord o1, AlignRecord o2)
    {
        return (0 == this.comp.compare(o1.record, o2.record)) ? true : false;
    }
}
