#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: GenomeAnalysis
   :synopsis: GenomeAnalysis implements genome analysis specifics
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
from libs.BUSCO.src.busco.analysis.BuscoAnalysis import BuscoAnalysis
from libs.BUSCO.src.busco.analysis.Analysis import NucleotideAnalysis, BLASTAnalysis
from libs.BUSCO.src.busco.busco_tools.prodigal import ProdigalRunner
from libs.BUSCO.src.busco.busco_tools.metaeuk import MetaeukRunner
from libs.BUSCO.src.busco.busco_tools.augustus import (
    AugustusRunner,
    GFF2GBRunner,
    NewSpeciesRunner,
    ETrainingRunner,
    OptimizeAugustusRunner,
)
from libs.BUSCO.src.busco.busco_tools.base import NoRerunFile, NoGenesError
from libs.BUSCO.src.busco.BuscoLogger import BuscoLogger
from libs.BUSCO.src.busco.BuscoLogger import LogDecorator as log
from abc import ABCMeta, abstractmethod
from configparser import NoOptionError
import time
import os
import glob
import pandas as pd
from collections import defaultdict
import subprocess
from libs.BUSCO.src.busco.Exceptions import BuscoError

logger = BuscoLogger.get_logger(__name__)


class GenomeAnalysis(NucleotideAnalysis, BuscoAnalysis, metaclass=ABCMeta):

    _mode = "genome"

    def __init__(self):
        super().__init__()
        #self.bbtools_runner = None

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()
        #if "genome" in self.config.get("busco_run", "mode"):
        #    # This avoids the child class euk_tran mode running bbtools
        #    self._run_bbtools()

    def init_tools(self):
        """
        Initialize tools needed for Genome Analysis.
        :return:
        """
        super().init_tools()
        #if "genome" in self.config.get("busco_run", "mode"):
        #    # This avoids the child class euk_tran mode running bbtools
        #    self.bbtools_runner = BBToolsRunner()

    #def _run_bbtools(self):
    #    self.bbtools_runner.configure_runner()
    #    if self.bbtools_runner.check_previous_completed_run():
    #        logger.info("Skipping BBTools run as already run")
    #    else:
    #        self.bbtools_runner.run()
    #    if len(self.bbtools_runner.metrics) == 0:
    #        self.bbtools_runner.parse_output()

    def reset(self):
        super().reset()
        #if self.bbtools_runner:
        #   self.bbtools_runner.reset()


class GenomeAnalysisProkaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a genome.
    """

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()
        self.prodigal_runner = None

    def cleanup(self):
        super().cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self._run_prodigal()
        self.run_hmmer(self.prodigal_runner.output_faa)
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        return

    def init_tools(self):
        """
        Init the tools needed for the analysis
        """
        super().init_tools()
        self.prodigal_runner = ProdigalRunner()

    def reset(self):
        super().reset()
        if (
            self.prodigal_runner
        ):  # If final run has already been run, then the prodigal_runner object in the final runner object will
            # still be set to None
            self.prodigal_runner.reset()

    @log("***** Run Prodigal on input to predict and extract genes *****", logger)
    def _run_prodigal(self):
        """
        Run Prodigal on input file to detect genes.
        :return:
        """
        self.prodigal_runner.configure_runner()
        if self.restart and self.prodigal_runner.check_previous_completed_run():
            logger.info("Skipping Prodigal run as it has already completed")
            self.prodigal_runner.run(restart=self.restart)
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.prodigal_runner.run()
        self.gene_details = self.prodigal_runner.gene_details
        self.sequences_nt = self.prodigal_runner.sequences_nt
        self.sequences_aa = self.prodigal_runner.sequences_aa

        return


class GenomeAnalysisEukaryotes(GenomeAnalysis):
    """
    This class runs a BUSCO analysis on a eukaryote genome.
    """

    def __init__(self):
        super().__init__()

        self.sequences_nt = {}
        self.sequences_aa = {}

    def cleanup(self):
        """
        This function cleans temporary files
        """
        super().cleanup()

    def init_tools(self):
        """
        Initialize all required tools for Genome Eukaryote Analysis:
        metaeuk
        :return:
        """
        super().init_tools()

        return

    def reset(self):
        super().reset()

    @abstractmethod
    def run_analysis(self):
        super().run_analysis()


class GenomeAnalysisEukaryotesAugustus(BLASTAnalysis, GenomeAnalysisEukaryotes):
    def __init__(self):
        super().__init__()
        self._long = self.config.getboolean("busco_run", "long")
        try:
            self._target_species_initial = self.config.get(
                "busco_run", "augustus_species"
            )
            self._target_species = self._target_species_initial
        except KeyError:
            raise BuscoError(
                "Something went wrong. Eukaryota datasets should specify an augustus species."
            )
        try:
            self._augustus_parameters = self.config.get(
                "busco_run", "augustus_parameters"
            ).replace(",", " ")
        except NoOptionError:
            self._augustus_parameters = ""
        self.mkblast_runner = None
        self.tblastn_runner = None
        self.augustus_runner = None
        self.gff2gb_runner = None
        self.new_species_runner = None
        self.etraining_runner = None
        self.optimize_augustus_runner = None

    def init_tools(self):
        super().init_tools()
        self.augustus_runner = AugustusRunner()
        self.gff2gb_runner = GFF2GBRunner()
        self.new_species_runner = NewSpeciesRunner()
        self.etraining_runner = ETrainingRunner()

        if self._long:
            self.optimize_augustus_runner = OptimizeAugustusRunner()

    def reset(self):
        super().reset()
        self.augustus_runner.reset()
        self.gff2gb_runner.reset()
        self.new_species_runner.reset()
        self.etraining_runner.reset()
        if self._long:
            self.optimize_augustus_runner.reset()

    def cleanup(self):
        try:
            self.augustus_runner.move_retraining_parameters()
            self.config.set(
                "busco_run", "augustus_species", self._target_species_initial
            )  # Reset parameter for batch mode
        except (OSError, AttributeError):
            pass
        super().cleanup()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        self._run_mkblast()
        self._run_tblastn()
        self._run_augustus(self.tblastn_runner.coords)
        self.gene_details = self.augustus_runner.gene_details
        self.run_hmmer(self.augustus_runner.output_sequences)
        self._rerun_analysis()

    def _rerun_augustus(self, coords):
        missing_and_fragmented_buscos = self.hmmer_runner.missing_buscos + list(
            self.hmmer_runner.fragmented_buscos.keys()
        )
        logger.info(
            "Re-running Augustus with the new metaparameters, number of target BUSCOs: {}".format(
                len(missing_and_fragmented_buscos)
            )
        )
        missing_and_fragmented_coords = {
            busco: coords[busco]
            for busco in coords
            if busco in missing_and_fragmented_buscos
        }
        logger.debug("Trained species folder is {}".format(self._target_species))
        self._run_augustus(missing_and_fragmented_coords)
        return

    @log(
        "Starting second step of analysis. The gene predictor Augustus is retrained using the results from the "
        "initial run to yield more accurate results.",
        logger,
    )
    def _rerun_analysis(self):

        self.augustus_runner.make_gff_files(self.hmmer_runner.single_copy_buscos)
        self._run_tblastn(
            missing_and_frag_only=True, ancestral_variants=self._has_variants_file
        )
        self._run_gff2gb()
        self._run_new_species()
        self.config.set(
            "busco_run", "augustus_species", self.new_species_runner.new_species_name
        )
        self._target_species = self.new_species_runner.new_species_name
        self._run_etraining()

        if self._long:
            self._run_optimize_augustus(self.new_species_runner.new_species_name)
            self._run_etraining()

        try:
            self._rerun_augustus(self.tblastn_runner.coords)
            self.gene_details.update(self.augustus_runner.gene_details)
            self.run_hmmer(self.augustus_runner.output_sequences)
            self.augustus_runner.make_gff_files(self.hmmer_runner.single_copy_buscos)
            self.augustus_runner.make_gff_files(self.hmmer_runner.multi_copy_buscos)
            self.augustus_runner.make_gff_files(self.hmmer_runner.fragmented_buscos)
            self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        except NoGenesError:
            logger.warning("No genes found on Augustus rerun.")

        return

    @log("Running Augustus gene predictor on BLAST search results.", logger)
    def _run_augustus(self, coords):
        self.augustus_runner.configure_runner(
            self.tblastn_runner.output_seqs,
            coords,
            self.sequences_aa,
            self.sequences_nt,
        )

        if self.restart and self.augustus_runner.check_previous_completed_run():
            run = "2nd" if self.augustus_runner.run_number == 2 else "1st"
            logger.info(
                "Skipping {} augustus run as output already processed".format(run)
            )
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.augustus_runner.run()
        self.augustus_runner.process_output()
        self.sequences_nt = self.augustus_runner.sequences_nt
        self.sequences_aa = self.augustus_runner.sequences_aa

    def _run_etraining(self):
        """Train on new training set (complete single copy buscos)"""
        self.etraining_runner.configure_runner(self.new_species_runner.new_species_name)
        if self.restart and self.etraining_runner.check_previous_completed_run():
            logger.info("Skipping etraining as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.etraining_runner.run()
        return

    @log("Converting predicted genes to short genbank files", logger)
    def _run_gff2gb(self):
        self.gff2gb_runner.configure_runner(self.hmmer_runner.single_copy_buscos)
        if self.restart and self.gff2gb_runner.check_previous_completed_run():
            logger.info("Skipping gff2gb conversion as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.gff2gb_runner.run()
        return

    @log(
        "All files converted to short genbank files, now training Augustus using Single-Copy Complete BUSCOs",
        logger,
    )
    def _run_new_species(self):
        """Create new species config file from template"""
        self.new_species_runner.configure_runner()
        if self.restart and self.new_species_runner.check_previous_completed_run():
            logger.info("Skipping new species creation as it has already been done")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.new_species_runner.run()
        return

    def _run_optimize_augustus(self, new_species_name):
        """ long mode (--long) option - runs all the Augustus optimization scripts (adds ~1 day of runtime)"""
        logger.warning(
            "Optimizing augustus metaparameters, this may take a very long time, started at {}".format(
                time.strftime("%m/%d/%Y %H:%M:%S")
            )
        )
        self.optimize_augustus_runner.configure_runner(
            self.augustus_runner.output_folder, new_species_name
        )
        self.optimize_augustus_runner.run()
        return


class GenomeAnalysisEukaryotesMetaeuk(GenomeAnalysisEukaryotes):
    def __init__(self):
        super().__init__()
        self.metaeuk_runner = None
        self.gene_details = {}
        self.gene_update_mapping = {}

    def init_tools(self):
        super().init_tools()

        self.metaeuk_runner = MetaeukRunner()

    def reset(self):
        super().reset()
        self.metaeuk_runner.reset()

    def run_analysis(self):
        """This function calls all needed steps for running the analysis."""
        super().run_analysis()
        incomplete_buscos = None
        for i in range(2):
            try:
                self._run_metaeuk(incomplete_buscos)
                self.gene_details.update(self.metaeuk_runner.gene_details)
                self.sequences_aa.update(self.metaeuk_runner.sequences_aa)
                self.sequences_nt.update(self.metaeuk_runner.sequences_nt)
                out = self.run_hmmer(
                    self.metaeuk_runner.pred_protein_seqs_modified,
                    busco_ids=incomplete_buscos,
                )
                incomplete_buscos = self.hmmer_runner.missing_buscos + list(
                    self.hmmer_runner.fragmented_buscos.keys()
                )
                if len(incomplete_buscos) == 0:
                    break
            except NoRerunFile:
                if i == 1:
                    logger.info("Metaeuk rerun did not find any genes")
                else:
                    raise NoGenesError("Metaeuk")

        try:
            self.metaeuk_runner.combine_run_results()
        except FileNotFoundError:
            # This exception should only happen if the rerun file does not exist. If the initial run file was
            # missing there would have been a BatchFatalError call above. The index 0 sets the "combined" file to the
            # output of the initial run.
            self.metaeuk_runner.combined_pred_protein_seqs = (
                self.metaeuk_runner.pred_protein_mod_files[0]
            )
            self.metaeuk_runner.combined_nucleotides_seqs = (
                self.metaeuk_runner.codon_mod_files[0]
            )
        self.hmmer_runner.write_buscos_to_file(self.sequences_aa, self.sequences_nt)
        self.write_gff_files()
        return out

    def write_gff_files(self):
        busco_seqs_folders = [
            self.hmmer_runner.single_copy_sequences_folder,
            self.hmmer_runner.multi_copy_sequences_folder,
            self.hmmer_runner.fragmented_sequences_folder,
        ]
        for folder in busco_seqs_folders:
            existing_gff_files = glob.glob(os.path.join(folder, "*.gff"))
            for f in existing_gff_files:
                os.remove(f)

        for gff_file in self.metaeuk_runner.gff_files[
            ::-1
        ]:  # Write GFF results to file, starting with the rerun results and then using the initial results
            self.metaeuk_runner.gff_file = gff_file
            self.metaeuk_runner.write_gff_files(*busco_seqs_folders)

    def _run_metaeuk(self, incomplete_buscos):
        self.metaeuk_runner.configure_runner(incomplete_buscos)
        if self.restart and self.metaeuk_runner.check_previous_completed_run():
            logger.info("Skipping Metaeuk run as already run")
        else:
            self.restart = False
            self.config.set("busco_run", "restart", str(self.restart))
            self.metaeuk_runner.run()

        self.metaeuk_runner.edit_result_headers()
        self.metaeuk_runner.get_gene_details()  # The gene details contain the overlaps that were removed when editing
        # the protein file, but it doesn't matter, as it is just a look-up
        # dictionary

    def validate_output(self):
        """Need to run this for both initial and rerun because it is possible that metaeuk matches overlap"""

        hmmer_results = self.hmmer_runner.merge_dicts()

        if len(hmmer_results) > 0:
            exon_records = self.get_exon_records(
                hmmer_results, self.hmmer_runner.run_number
            )
            df = self.exons_to_df(exon_records)
            overlaps = self.find_overlaps(df)
            if len(overlaps) > 0:
                discarded_exon_lengths = self.handle_overlaps(overlaps, df)
                # df.drop(inds_to_remove, inplace=True)
                complete, matched_genes_complete = self.reconstruct_hmmer_results(
                    df, discarded_exon_lengths, self.hmmer_runner.is_complete
                )
                v_large, matched_genes_v_large = self.reconstruct_hmmer_results(
                    df, discarded_exon_lengths, self.hmmer_runner.is_very_large
                )
                fragmented, matched_genes_fragmented = self.reconstruct_hmmer_results(
                    df, discarded_exon_lengths, self.hmmer_runner.is_fragment
                )

                # Update hmmer runner with new dictionaries
                self.hmmer_runner.is_complete = complete
                self.hmmer_runner.is_very_large = v_large
                self.hmmer_runner.is_fragment = fragmented
                self.hmmer_runner.matched_genes_complete = matched_genes_complete
                self.hmmer_runner.matched_genes_vlarge = matched_genes_v_large
                self.hmmer_runner.matched_genes_fragment = matched_genes_fragmented
                self.hmmer_runner.gene_details = self.gene_details
        return

    def get_exon_records(
        self, busco_dict, run_number
    ):  # Placed in the GenomeAnalysis module because it draws on both hmmer_runner and metaeuk_runner methods

        initial_run_results = self.metaeuk_runner.headers_files[0]
        if run_number == 2:
            rerun_results = self.metaeuk_runner.headers_files[1]

        exon_records = []
        for busco_id, gene_match in busco_dict.items():
            for gene_id, details in gene_match.items():
                sequence, coords = details[0]["orig gene ID"].rsplit(":", 1)
                gene_start, gene_end = coords.split("-")
                strand = self.gene_details[gene_id][0]["strand"]
                score = details[0]["bitscore"]

                # Need to determine run using HMMER results instead of metaeuk results. This is because exons can be
                # matched identically to different BUSCO IDs, both on the same and on different runs. The presence of a
                # match in the metaeuk rerun results does not indicate that the HMMER match in question is associated
                # with that metaeuk match
                run_found = (
                    "2"
                    if run_number == 2
                    and os.path.exists(
                        os.path.join(
                            self.hmmer_runner.rerun_results_dir,
                            "{}.out".format(busco_id),
                        )
                    )
                    else "1"
                )
                try:
                    if run_found == "2":
                        matches = subprocess.check_output(
                            [
                                "grep",
                                "{}|{}|.*|{}|{}|".format(
                                    sequence, strand, gene_start, gene_end
                                ),
                                rerun_results,
                            ]
                        ).decode("utf-8")
                    else:
                        matches = subprocess.check_output(
                            [
                                "grep",
                                "{}|{}|.*|{}|{}|".format(
                                    sequence, strand, gene_start, gene_end
                                ),
                                initial_run_results,
                            ]
                        ).decode("utf-8")
                except subprocess.CalledProcessError:
                    raise BuscoError(
                        "Unable to parse metaeuk results. This typically occurs because sequence headers "
                        "contain pipes ('|'). Metaeuk uses pipes as delimiters in the results files. The "
                        "additional pipes interfere with BUSCO's ability to accurately parse the results."
                        "To fix this problem remove any pipes from sequence headers and try again."
                    )

                # The following line is a relic from when the previous grep search tried to match busco_id instead of
                # gene details. The find_match method is still needed though to clean up the match, even though it
                # redundantly matches the gene coordinates again.
                good_match = self.metaeuk_runner.find_match(
                    matches,
                    ["|{}|".format(gene_start), "|{}|".format(gene_end), sequence],
                )

                if good_match:
                    low_coords, high_coords = self.metaeuk_runner.extract_exon_coords(
                        good_match
                    )
                    if low_coords[0] > high_coords[0]:  # for negative strand exons the order is reversed
                        low_coords, high_coords = high_coords, low_coords
                    trimmed_low = int(gene_id.split(":")[-1].split("-")[0])
                    trimmed_high = int(gene_id.split(":")[-1].split("-")[1])
                    for i, entry in enumerate(low_coords):
                        if int(entry) < trimmed_low or int(entry) > trimmed_high:  # don't include exons that were previously removed due to overlaps
                            continue
                        record = (
                            busco_id,
                            sequence,
                            entry,
                            high_coords[i],
                            strand,
                            score,
                            run_found,
                            details[0]["orig gene ID"],
                        )
                        exon_records.append(record)
        return exon_records

    def reconstruct_hmmer_results(self, df, discarded_exon_lengths, hmmer_result_dict):
        busco_groups = df.groupby(["Busco id"])
        inds_to_remove = set(discarded_exon_lengths.keys())
        hmmer_result_dict_new = defaultdict(dict)
        matched_genes_new = defaultdict(list)
        for busco_id, matches in hmmer_result_dict.items():
            try:
                busco_group = busco_groups.get_group(busco_id)
            except KeyError:  # if busco was removed during overlap filtering
                continue
            busco_gene_groups = busco_group.groupby("Orig gene ID")
            for gene_match, busco_gene_group in busco_gene_groups:
                if gene_match not in matches:
                    if gene_match not in self.gene_update_mapping:
                        continue
                    else:
                        gene_match = self.gene_update_mapping[gene_match]
                new_gene_start = gene_match.split(":")[-1].split("-")[
                    0
                ]  # these two lines are not really used - they just initialize values that will be changed
                new_gene_stop = gene_match.split(":")[-1].split("-")[1]
                start_trim = 0
                end_trim = 0
                group_indices = set(busco_gene_group.index)
                intersection = group_indices.intersection(inds_to_remove)
                if len(intersection) > 0:
                    if intersection == group_indices:
                        continue  # remove entire gene - don't add to new dict
                    ordered_exons = busco_gene_group.sort_values(by="Start").reset_index()
                    new_indices = ordered_exons.index
                    seq = ordered_exons.loc[0]["Sequence"]

                    for idx in new_indices:
                        old_index = ordered_exons.loc[idx]["index"]
                        if old_index in intersection:
                            start_trim += discarded_exon_lengths[old_index]
                        else:
                            new_gene_start = df.loc[old_index]["Start"]
                            break
                    for idx in new_indices[::-1]:
                        old_index = ordered_exons.loc[idx]["index"]
                        if old_index in intersection:
                            end_trim += discarded_exon_lengths[old_index]
                        else:
                            new_gene_stop = df.loc[old_index]["Stop"]
                            break
                    new_gene_match = "{}:{}-{}".format(
                        seq, new_gene_start, new_gene_stop
                    )
                else:
                    new_gene_match = gene_match

                details = matches[gene_match]
                df_strand = busco_gene_group["Strand"].iloc[0]
                hmmer_result_dict_new[busco_id].update({new_gene_match: details})
                matched_genes_new[new_gene_match].append(busco_id)
                self.gene_details[new_gene_match] = [
                    {
                        "gene_start": new_gene_start,
                        "gene_end": new_gene_stop,
                        "strand": df_strand,
                    }
                ]
                if new_gene_match != gene_match:
                    trimmed_sequence_aa, trimmed_sequence_nt = self.trim_sequence(
                        gene_match, start_trim, end_trim
                    )
                    self.gene_update_mapping[gene_match] = new_gene_match
                else:
                    try:
                        trimmed_sequence_aa = self.metaeuk_runner.sequences_aa[gene_match]
                        trimmed_sequence_nt = self.metaeuk_runner.sequences_nt[gene_match]
                    except KeyError:  # happens on the second run if the first run trimmed the sequence already
                        trimmed_sequence_aa = self.sequences_aa[gene_match]
                        trimmed_sequence_nt = self.sequences_nt[gene_match]
                self.sequences_aa[new_gene_match] = trimmed_sequence_aa
                self.sequences_nt[new_gene_match] = trimmed_sequence_nt
        return hmmer_result_dict_new, matched_genes_new

    def trim_sequence(self, old_gene_match, start_trim, end_trim):
        old_sequence_aa = self.metaeuk_runner.sequences_aa[old_gene_match]
        old_sequence_nt = self.metaeuk_runner.sequences_nt[old_gene_match]

        new_sequence_nt = old_sequence_nt[start_trim : len(old_sequence_nt) - end_trim]
        if start_trim % 3 != 0 and end_trim % 3 != 0:
            raise BuscoError(
                "Problem reconstructing amino acid sequence after extracting exons"
            )
        new_sequence_aa = old_sequence_aa[
            int(start_trim / 3) : len(old_sequence_aa) - int(end_trim / 3)
        ]
        return new_sequence_aa, new_sequence_nt

    def exons_to_df(self, records):
        if self._mode == "genome":
            logger.info("Validating exons and removing overlapping matches")

        df = self.metaeuk_runner.records_to_df(records)
        df["Start"] = df["Start"].astype(int)
        df["Stop"] = df["Stop"].astype(int)
        df["Score"] = df["Score"].astype(float)
        df["Run found"] = df["Run found"].astype(int)
        return df

    def find_overlaps(self, df):
        overlaps = self.metaeuk_runner.test_for_overlaps(df, sort=True)
        logger.info("{} candidate overlapping regions found".format(len(overlaps)))
        logger.info("{} exons in total".format(len(df)))
        return overlaps

    def handle_overlaps(self, overlaps, df):
        discarded_exon_lengths = {}
        # Consider three gene matches, A, B and C, with decreasing scores respectively.
        # If A overlaps B and B overlaps C but A does not overlap C then the order they are treated affects the outcome.
        for n, overlap_inds in enumerate(overlaps):
            if (
                overlap_inds[0] in discarded_exon_lengths
                or overlap_inds[1] in discarded_exon_lengths
            ):
                continue
            else:
                # logger.info("Overlap {}:".format(n))
                bad_inds = self.handle_diff_busco_overlap(overlap_inds, df)
                for idx in bad_inds:
                    discarded_exon_lengths[idx] = (
                        abs(df.iloc[idx]["Stop"] - df.iloc[idx]["Start"]) + 1
                    )
        return discarded_exon_lengths

    def handle_diff_busco_overlap(self, overlap_inds, df):
        match1 = df.loc[overlap_inds[0]]
        match2 = df.loc[overlap_inds[1]]
        seq = match1["Sequence"]
        busco_match1 = match1["Busco id"]
        gene_match1 = match1["Orig gene ID"]
        run_match1 = match1["Run found"]
        busco_match2 = match2["Busco id"]
        gene_match2 = match2["Orig gene ID"]
        run_match2 = match2["Run found"]
        exons1 = df.loc[
            (df["Busco id"] == busco_match1)
            & (df["Sequence"] == seq)
            & (df["Orig gene ID"] == gene_match1)
        ]
        exons2 = df.loc[
            (df["Busco id"] == busco_match2)
            & (df["Sequence"] == seq)
            & (df["Orig gene ID"] == gene_match2)
        ]
        hmmer_run_folder1 = (
            self.hmmer_runner.initial_results_dir
            if run_match1 == 1
            else self.hmmer_runner.rerun_results_dir
        )
        hmmer_run_folder2 = (
            self.hmmer_runner.initial_results_dir
            if run_match2 == 1
            else self.hmmer_runner.rerun_results_dir
        )
        hmmer_result1 = os.path.join(hmmer_run_folder1, "{}.out".format(busco_match1))
        hmmer_result2 = os.path.join(hmmer_run_folder2, "{}.out".format(busco_match2))
        hmmer_match_details1 = self.hmmer_runner.parse_hmmer_output(
            hmmer_result1, busco_match1
        )
        hmmer_match_details2 = self.hmmer_runner.parse_hmmer_output(
            hmmer_result2, busco_match2
        )
        gene_id1 = match1["Orig gene ID"]
        gene_id2 = match2["Orig gene ID"]
        if (
            hmmer_match_details1[gene_id1]["score"]
            > hmmer_match_details2[gene_id2]["score"]
        ):
            priority_match = hmmer_match_details1
            secondary_match = hmmer_match_details2
            priority_exons = exons1
            secondary_exons = exons2
            priority_gene_id = gene_id1
            secondary_gene_id = gene_id2
        else:
            priority_match = hmmer_match_details2
            secondary_match = hmmer_match_details1
            priority_exons = exons2
            secondary_exons = exons1
            priority_gene_id = gene_id2
            secondary_gene_id = gene_id1

        priority_env_coords = iter(priority_match[priority_gene_id]["env_coords"])
        secondary_env_coords = iter(secondary_match[secondary_gene_id]["env_coords"])
        priority_used_exons, priority_unused_exons = self.find_unused_exons(
            priority_env_coords, priority_exons
        )
        secondary_used_exons, secondary_unused_exons = self.find_unused_exons(
            secondary_env_coords, secondary_exons
        )

        priority_used_exons = (
            pd.DataFrame.from_records(priority_used_exons, index="index")
            if priority_used_exons
            else None
        )
        priority_unused_exons = (
            pd.DataFrame.from_records(priority_unused_exons, index="index")
            if priority_unused_exons
            else None
        )
        secondary_used_exons = (
            pd.DataFrame.from_records(secondary_used_exons, index="index")
            if secondary_used_exons
            else None
        )
        secondary_unused_exons = (
            pd.DataFrame.from_records(secondary_unused_exons, index="index")
            if secondary_unused_exons
            else None
        )

        indices_to_remove = []
        # Check if secondary match uses priority match exons
        used_exons = pd.concat([priority_used_exons, secondary_used_exons])
        overlaps = self.metaeuk_runner.test_for_overlaps(used_exons)

        if len(overlaps) > 0:
            # Remove secondary match
            indices_to_remove = secondary_exons.index
            return indices_to_remove

        # Check to see if unused priority exons are used by secondary match
        indices_to_remove.extend(
            self.get_indices_to_remove(secondary_used_exons, priority_unused_exons)
        )

        # Check to see if unused secondary exons are used by priority match
        indices_to_remove.extend(
            self.get_indices_to_remove(priority_used_exons, secondary_unused_exons)
        )

        # Check to see if unused secondary exons overlap with priority unused exons
        indices_to_remove.extend(
            self.get_indices_to_remove(priority_unused_exons, secondary_unused_exons)
        )

        return indices_to_remove

    def get_indices_to_remove(self, priority_exons, secondary_exons):
        indices_to_remove = []
        try:
            exons_to_check = pd.concat([priority_exons, secondary_exons])
        except ValueError:
            # all exons are None
            return indices_to_remove

        overlaps = self.metaeuk_runner.test_for_overlaps(exons_to_check)
        if len(overlaps) > 0:
            for overlap in overlaps:
                if overlap[0] in indices_to_remove or overlap[1] in indices_to_remove:
                    continue
                match1 = exons_to_check.loc[overlap[0]]
                if secondary_exons is not None:
                    index_to_remove = (
                        overlap[0]
                        if secondary_exons.iloc[0]["Busco id"] == match1["Busco id"]
                        else overlap[1]
                    )  # It is possible that secondary_exons is None so need to check before using iloc
                else:
                    match2 = exons_to_check.loc[overlap[1]]
                    index_to_remove = (
                        overlap[0] if match1["Score"] < match2["Score"] else overlap[1]
                    )

                exons_to_remove = secondary_exons if index_to_remove in secondary_exons.index else priority_exons
                indices_to_remove.extend(list(exons_to_remove.index))
        return indices_to_remove

    def find_unused_exons(self, env_coords, exons):
        remaining_hmm_region = 0
        unused_exons = []
        used_exons = []
        hmm_coords = next(env_coords)
        exon_cumul_len = 0
        for idx, entry in exons.iterrows():
            entry["index"] = idx
            exon_matched = False
            exon_size_nt = int(entry["Stop"]) - int(entry["Start"]) + 1
            if not exon_size_nt % 3 == 0:
                raise BuscoError(
                    "The exon coordinates contain fractional reading frames and are ambiguous."
                )
            exon_size_aa = exon_size_nt / 3
            exon_cumul_len += exon_size_aa
            if remaining_hmm_region > exon_size_aa:
                remaining_hmm_region -= exon_size_aa
                exon_matched = True

            elif remaining_hmm_region:
                exon_matched = True

            elif hmm_coords:
                while hmm_coords[0] < exon_cumul_len + 1:
                    # hmm starts within exon
                    exon_matched = True
                    if hmm_coords[1] <= exon_cumul_len + 1:
                        # hmm ends within exon; cycle to the next hmm region
                        try:
                            hmm_coords = next(env_coords)
                        except StopIteration:
                            hmm_coords = None
                            break
                        continue
                    else:
                        remaining_hmm_region = hmm_coords[1] - exon_size_aa + 1
                        break
            if exon_matched:
                used_exons.append(entry)
            else:
                unused_exons.append(entry)
        used_exons, unused_exons = self.adjust_exon_categories(used_exons, unused_exons)
        return used_exons, unused_exons

    @staticmethod
    def adjust_exon_categories(used_exons, unused_exons):
        """
        Ensure there are no unused exons sandwiched between used exons
        :param used_exons:
        :param unused_exons:
        :return:
        """

        used_exons_start = [x["Start"] for x in used_exons]
        used_exons_end = [x["Stop"] for x in used_exons]
        start = min(used_exons_start)
        stop = max(used_exons_end)
        exons_to_remove = set()
        unused_indices = [exon["index"] for exon in unused_exons]
        for exon in unused_exons:
            if not exon["index"] in exons_to_remove and (
                (exon["Start"] >= start and exon["Start"] < stop)
                or (exon["Stop"] > start and exon["Stop"] < stop)
            ):
                # find exons that either start or stop within the "used" range
                used_exons.append(exon)
                exons_to_remove.add(exon["index"])
        for idx in list(exons_to_remove):
            idx2 = unused_indices.index(idx)
            unused_exons.pop(idx2)
            unused_indices.pop(idx2)

        return used_exons, unused_exons

    def cleanup(self):
        try:
            self.metaeuk_runner.remove_tmp_files()
        except OSError:
            pass
        super().cleanup()
