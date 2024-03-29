# Original code written by: Thomas Murray and taken from https://github.com/tom-murray98/minFQ_BAM/blob/master/BAM_RG_tag_extractor.py

import pysam


class ReadBam:
    """
    :param bam_file: A bam format alignment file.
    :return:
        id_tag : ID
        dt_tag: time of run
        basecall_model_tag: Basecall Model ID
        runid_tag: run ID
        lb: sample ID
        pl: outputs 'ONT', device organisation?
        pm: device position
        pu: flow cell ID
        al: unknown (outputs unclassified when present)
    """

    def __init__(self, bam_file=None):
        self.sam_file = None
        self.bam_file = bam_file
        self.mapped_reads = 0
        self.unmapped_reads = 0
        if "pass" in self.bam_file:
            self.state = "pass"
        else:
            self.state = "fail"


    def get_rg_tags(self):
        rg_tags = self.sam_file.header.get("RG", [])
        if not rg_tags:
            print("This BAM file does not contain an @RG field")
            return None
        else:
            for rg_tag in rg_tags:
                id_tag = rg_tag.get("ID", None)
                dt_tag = rg_tag.get("DT", None)
                ds_tag = rg_tag.get("DS", None)
                ds_tags = ds_tag.split(" ")
                basecall_model_tag = ds_tags[1].replace("basecall_model=", "")
                runid_tag = ds_tags[0].replace("runid=", "")
                lb_tag = rg_tag.get("LB", None)
                pl_tag = rg_tag.get("PL", None)
                pm_tag = rg_tag.get("PM", None)
                pu_tag = rg_tag.get("PU", None)
                al_tag = rg_tag.get("al", None)

                return (
                    id_tag,
                    dt_tag,
                    basecall_model_tag,
                    runid_tag,
                    lb_tag,
                    pl_tag,
                    pm_tag,
                    pu_tag,
                    al_tag,
                )

    def read_bam(self):
        if self.bam_file:
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)

            (
                id_tag,
                dt_tag,
                basecall_model_tag,
                runid_tag,
                lb_tag,
                pl_tag,
                pm_tag,
                pu_tag,
                al_tag,
            ) = self.get_rg_tags()

            yield (
                # name,
                # seq,
                # qual,
                # start_time_pr,
                # channel,
                # read_number,
                # read_basecall_id,
                id_tag,
                dt_tag,
                basecall_model_tag,
                runid_tag,
                lb_tag,
                pl_tag,
                pm_tag,
                pu_tag,
                al_tag,
            )

    def process_reads(self):
        for (
            # name,
            # seq,
            # qual,
            # start_time_pr,
            # channel,
            # read_number,
            # read_basecall_id,
            id_tag,
            dt_tag,
            basecall_model_tag,
            runid_tag,
            lb_tag,
            pl_tag,
            pm_tag,
            pu_tag,
            al_tag,
        ) in self.read_bam():
            bam_read = {
                "ID": id_tag,
                "time_of_run": dt_tag,
                "sample_id": lb_tag,
                "basecall_model": basecall_model_tag,
                "runid": runid_tag,
                "platform": pl_tag,
                "flow_cell_id": pu_tag,
                "device_position": pm_tag,
                "al": al_tag,
                "state": self.state,
                # "read_id": name,
                # "sequence": seq,
                # "sequence_length": len(str(seq)),
                # "quality": qual,
                # "start_time_per_run": start_time_pr,
                # "channel": channel,
                # "read_number": read_number,
                # "read_basecall_id": read_basecall_id,
            }

            filtered_reads_count = 0
            yield_tracking = 0
            readset = set()
            for read in self.sam_file.fetch():
                # Check if read is not unmapped and not a secondary alignment
                if not read.is_unmapped and not read.is_secondary:
                    filtered_reads_count += 1
                    if read.query_name not in readset:
                        readset.add(read.query_name)
                        #ToDo: Should this only happen for each read once? Or should it happen every time an alignment is seen?
                        #ToDo: Strictly this is aligned length and not total length??
                        if read.infer_query_length() and read.infer_query_length()>0:
                            yield_tracking += read.infer_query_length()
                    #print(f"Read name: {read.query_name}")
            self.mapped_reads=len(readset)
            self.yield_tracking = yield_tracking
            self.unmapped_reads=self.sam_file.unmapped
            #print(f"Mapped reads:{filtered_reads_count}")
            #print(f"Total reads:{filtered_reads_count + self.sam_file.unmapped}")
            return bam_read
