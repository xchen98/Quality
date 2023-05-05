#!/usr/bin/env python3
# coding: utf-8

"""
.. module:: BuscoPlacer
   :synopsis: BuscoPlacer implements methods required for automatically selecting the appropriate dataset
   to be used during BUSCO analysis
.. versionadded:: 4.0.0
.. versionchanged:: 5.4.0

Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""
import json
import os
import re
from libs.BUSCO.src.busco.BuscoLogger import BuscoLogger
from libs.BUSCO.src.busco.BuscoLogger import LogDecorator as log
from libs.BUSCO.src.busco.Exceptions import BuscoError
from Bio import SeqIO
from libs.BUSCO.src.busco.busco_tools.sepp import SEPPRunner

logger = BuscoLogger.get_logger(__name__)


class NoMarkersError(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        return self.value


class BuscoPlacer:

    _logger = BuscoLogger.get_logger(__name__)

    @log(
        "***** Searching tree for chosen lineage to find best taxonomic match *****\n",
        logger,
    )
    def __init__(self, config, run_folder, protein_seqs, single_copy_buscos):
        self._config = config
        self._params = config
        self.mode = self._config.get("busco_run", "mode")
        self.cpus = self._config.get("busco_run", "cpu")
        self.restart = self._config.getboolean("busco_run", "restart")
        self.run_folder = run_folder
        self.placement_folder = os.path.join(run_folder, "placement_files")
        if self.restart:
            os.makedirs(self.placement_folder, exist_ok=True)
        else:
            os.mkdir(self.placement_folder)
        self.downloader = self._config.downloader
        self.datasets_version = self._config.get("busco_run", "datasets_version")
        self.protein_seqs = protein_seqs
        self.single_copy_buscos = single_copy_buscos  # dict
        self.init_tools()

    def _download_placement_files(self):
        self.ref_markers_file = self.downloader.get(
            "list_of_reference_markers.{0}_{1}.txt".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        self.tree_nwk_file = self.downloader.get(
            "tree.{0}_{1}.nwk".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        self.tree_metadata_file = self.downloader.get(
            "tree_metadata.{0}_{1}.txt".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        self.supermatrix_file = self.downloader.get(
            "supermatrix.aln.{0}_{1}.faa".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        self.taxid_busco_file = self.downloader.get(
            "mapping_taxids-busco_dataset_name.{0}_{1}.txt".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        self.taxid_lineage_file = self.downloader.get(
            "mapping_taxid-lineage.{0}_{1}.txt".format(
                os.path.basename(self.run_folder).split("_")[-2], self.datasets_version
            ),
            "placement_files",
        )
        return

    def _get_placement_file_versions(self):
        placement_file_versions = [
            os.path.basename(filepath)
            for filepath in [
                self.ref_markers_file,
                self.tree_nwk_file,
                self.tree_metadata_file,
                self.supermatrix_file,
                self.taxid_busco_file,
                self.taxid_lineage_file,
            ]
        ]
        return placement_file_versions

    @log("Extract markers...", logger)
    def define_dataset(self):
        # If mode is genome, substitute input with prodigal/augustus output
        self._download_placement_files()
        placement_file_versions = self._get_placement_file_versions()
        try:
            self._extract_marker_sequences()
            self._run_sepp()

            dataset = self._pick_dataset()
        except NoMarkersError:
            root_lineage = self._config.get("busco_run", "name")
            logger.info(
                "No marker genes were found. Root lineage {} is kept".format(
                    root_lineage
                )
            )
            dataset = (root_lineage, None, None)

        return dataset, placement_file_versions

    def init_tools(self):
        setattr(SEPPRunner, "config", self._config)
        self.sepp_runner = SEPPRunner()

    def _pick_dataset(self):

        run_folder = self.run_folder

        # load busco dataset name by id in a dict {taxid:name}
        datasets_mapping = {}

        with open(self.taxid_busco_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                tax_id = parts[0]
                dataset = parts[1].split(",")[0]
                datasets_mapping.update({tax_id: dataset})

        # load the lineage for each taxid in a dict {taxid:reversed_lineage}
        # lineage is 1:2:3:4:5:6 => {6:[6,5,4,3,2,1]}
        lineages = set()
        parents = {}
        taxid_dataset = {}
        for t in datasets_mapping:
            taxid_dataset.update({t: t})
        with open(self.taxid_lineage_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                lineage = line.strip().split("\t")[4]
                lineages.add(lineage)

                # for each line, e.g. 6\t1:2:3:4:5:6, create/update the lineage for each level
                # 6:[1,2,3,4,5,6], 5:[1,2,3,4,5], 4:[1,2,3,4], etc.
                levels = lineage.split(",")

                for i, t in enumerate(levels):
                    parents.update({t: levels[0 : i + 1][::-1]})

        for t in parents:
            for p in parents[t]:  # get the deepest parent, not the root one
                if p in datasets_mapping:
                    taxid_dataset.update({t: p})
                    break
        # load json
        # load "tree" in a string
        # load placements
        # obtain a dict of taxid num of markers
        # figure out which taxid to use by using the highest number of markers and some extra rules

        try:
            with open(
                os.path.join(self.placement_folder, "output_placement.json")
            ) as json_file:
                data = json.load(json_file)
            tree = data["tree"]
            placements = data["placements"]
        except FileNotFoundError:
            raise BuscoError(
                "Placements failed. Try to rerun increasing the memory or select a lineage manually."
            )

        node_weight = {}
        n_p = 0
        for placement in placements:
            n_p += 1
            for individual_placement in placement["p"]:
                # find the taxid in tree
                node = individual_placement[0]

                match = re.findall(  # deal with weird character in the json file, see the output yourself.
                    # if this pattern is inconsistant with pplacer version, it may break buscoplacer.
                    "[^0-9][0-9]*:[0-9]*[^0-9]{0,1}[0-9]*[^0-9]{0,2}[0-9]*\[%s\]"
                    % node,
                    tree,
                )
                # extract taxid:
                try:
                    if re.match("^[A-Za-z]", match[0]):
                        taxid = match[0][7:].split(":")[0]
                    else:
                        taxid = match[0][1:].split(":")[0]
                except IndexError as e:
                    raise e
                if taxid_dataset[taxid] in node_weight:
                    node_weight[taxid_dataset[taxid]] += 1
                else:
                    node_weight[taxid_dataset[taxid]] = 1
                break  # Break here to keep only the best match. In my experience, keeping all does not change much.
        type(self)._logger.debug("Placements counts by node are: %s" % node_weight)

        # from here, define which placement can be trusted
        max_markers = 0
        choice = []

        for key in node_weight:
            type(self)._logger.debug(
                "%s markers assigned to the taxid %s" % (node_weight[key], key)
            )

        # taxid for which no threshold or minimal amount of placement should be considered.
        # If it is the best, go for it.
        no_rules = ["204428"]

        ratio = 2.5
        if run_folder.split("/")[-1].split("_")[-2] == "archaea":
            ratio = 1.2
        min_markers = 12

        node_with_max_markers = None
        for n in node_weight:
            if node_weight[n] > max_markers:
                max_markers = node_weight[n]
                node_with_max_markers = n
        if node_with_max_markers in no_rules:
            choice = [node_with_max_markers]
        else:
            for n in node_weight:
                # if the ration between the best and the current one is not enough, keep both
                if node_weight[n] * ratio >= max_markers:
                    choice.append(n)
        if len(choice) > 1:
            # more than one taxid should be considered, pick the common ancestor
            choice = self._get_common_ancestor(choice, parents)
        elif len(choice) == 0:
            if run_folder.split("/")[-1].split("_")[-2] == "bacteria":
                choice.append("2")
            elif run_folder.split("/")[-1].split("_")[-2] == "archaea":
                choice.append("2157")
            elif run_folder.split("/")[-1].split("_")[-2] == "eukaryota":
                choice.append("2759")
        if max_markers < min_markers and not (choice[0] in no_rules):
            if run_folder.split("/")[-1].split("_")[-2] == "bacteria":
                key_taxid = "2"
            elif run_folder.split("/")[-1].split("_")[-2] == "archaea":
                key_taxid = "2157"
            elif run_folder.split("/")[-1].split("_")[-2] == "eukaryota":
                key_taxid = "2759"
            else:
                key_taxid = None  # unexpected. Should throw an exception or use assert.
            type(self)._logger.info(
                "Not enough markers were placed on the tree (%s). Root lineage %s is kept"
                % (max_markers, datasets_mapping[taxid_dataset[key_taxid]])
            )
            lineage = datasets_mapping[taxid_dataset[key_taxid]]

        else:
            type(self)._logger.info(
                "Lineage %s is selected, supported by %s markers out of %s"
                % (
                    datasets_mapping[taxid_dataset[choice[0]]],
                    max_markers,
                    sum(node_weight.values()),
                )
            )
            lineage = datasets_mapping[taxid_dataset[choice[0]]]
        lineage = "{}_{}".format(
            lineage, self._config.get("busco_run", "datasets_version")
        )
        placed_markers = sum(node_weight.values())
        return [
            lineage,
            max_markers,
            placed_markers,
        ]

    @staticmethod
    def _get_common_ancestor(choice, parents):
        # starts with the parents of the first choice
        all_ancestors = set(parents[choice[0]])
        # order will be lost with sets, so keep in a list the lineage of one entry to later pick the deepest ancestor
        ordered_lineage = []
        for c in choice:
            if len(parents[c]) > len(ordered_lineage):
                # probably useless. Init with parents[choice[0] should work
                ordered_lineage = parents[c]
            # keep in set only entries that are in the currently explored lineage
            all_ancestors = all_ancestors.intersection(parents[c])

        # go through the ordered list of the deepest linage until you found a common ancestor.
        for parent in ordered_lineage:
            if parent in all_ancestors:
                return [parent]

    @log("Place the markers on the reference tree...", logger)
    def _run_sepp(self):
        self.sepp_runner.configure_runner(
            self.tree_nwk_file,
            self.tree_metadata_file,
            self.supermatrix_file,
            self.downloader,
        )
        if self.restart and self.sepp_runner.check_previous_completed_run():
            logger.info("Skipping SEPP run as it has already been completed")
        else:
            self.restart = False
            self._config.set("busco_run", "restart", str(self.restart))
            self.sepp_runner.run()
            self.sepp_runner.cleanup()

    def _extract_marker_sequences(self):
        """
        This function extracts all single copy BUSCO genes from a protein run folder
        """

        with open(self.ref_markers_file, "r") as f:
            marker_list = [line.strip() for line in f]

        marker_genes_names = []
        for busco, gene_matches in self.single_copy_buscos.items():
            if busco in marker_list:
                marker_genes_names.append(
                    list(gene_matches.keys())[0]
                )  # The list should only have one entry because they are single copy buscos

        if len(marker_genes_names) == 0:
            raise NoMarkersError

        marker_genes_records = []
        if isinstance(self.protein_seqs, (str,)):
            list_protein_seqs = [self.protein_seqs]
        else:
            list_protein_seqs = self.protein_seqs

        for protein_seqs in list_protein_seqs:
            with open(protein_seqs, "r") as prot_seqs:
                for record in SeqIO.parse(prot_seqs, "fasta"):
                    if record.id in marker_genes_names:
                        record.seq = record.seq.rstrip("*")
                        record.description = ""
                        marker_genes_records.append(record)

        marker_genes_file = os.path.join(self.placement_folder, "marker_genes.fasta")
        with open(marker_genes_file, "w") as output:
            SeqIO.write(marker_genes_records, output, "fasta")
