from abc import ABC, abstractmethod
import logging
from .core import Region, SignalSegment


class RegionParserBase(ABC):
    """Base class for region file parsers"""
    
    @abstractmethod
    def parse(self, fp):
        """Read and parse a file containing genomic regions

        Parameters
        ----------
        fp : str
            Path to the file.

        Returns
        -------
        List[Regions]
            List of regions in the file
        """
        pass


class BedParser(RegionParserBase):
    """Parser for bed files.

    The parser reads optional fields. If not present in the file defaults
    will be used:
    - name: None
    - score: 0
    - strand: +
    """
    
    def parse(self, fp):
        result = []
        log_missing_info = False
        
        with open(fp) as f:
            for line in f:
                line = line.split()
                if line[0].lower() == 'track':  # ignoring track lines
                    continue
                chromosome, start, end = line[:3]
                start = int(start)
                end = int(end)
                name = None
                score = 0
                pos_strand = True
                
                try:
                    name = line[3]
                    score = float(line[4])
                    pos_strand = line[5] == '+'
                except IndexError:
                    log_missing_info = True
                
                result.append(Region(chromosome, start, end, name, pos_strand,
                                     score))
        
        if log_missing_info:
            logging.info('Some optional fields were missing in the file')
        
        return result


def read_bedgraph(fp):
    """Read and parse a bedgraph file.
    
    Parameters
    ----------
    fp : str
        Path to the file

    Returns
    -------
    list SignalSegments
        Signal from the bedgraph file.
    """
    result = []
    with open(fp) as f:
        for line in f:
            line = line.split()
            if line[0] == 'track':
                continue
            chromosome, start, end, value = line
            result.append(SignalSegment(chromosome, int(start), int(end),
                                        float(value)))
    return result
