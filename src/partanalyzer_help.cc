/*
*    This file is part of Partanalyzer.
*
*    Partanalyzer is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Partanalyzer is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Partanalyzer.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/

/*
#include "partanalyzer.h"
#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"
*/

#ifndef _HELPFUNCTIONS
#define _HELPFUNCTIONS 1

#include "partanalyzer_help.h"

void printCopyright(){
        cout<<"# "<<program<<" Version "<<VERSION<<endl;
        cout<<"# Copyright (c) Miguel A. Santos, May. 2008-2009 .Build "<<__DATE__<<endl;
        cout<<"# Licensed under the GNU GPL version 3 or later.       "<<endl;
        cout<<"# (see http://www.gnu.org/copyleft/gpl.html )          "<<endl;
        cout<<"#"<<endl;
}
void printVersion(){
        cout<<"# "<<program<<" Version "<<VERSION<<". Build "<<__DATE__<<endl;
}
void printCommandLineError(const string label){ //Defaults: label=""
        cout<<_programb_<<" : "<<label<<" : Unknown option or incorrect syntax. Use -h for syntax summary."<<endl;
	exit(1);
}

void printCommandLineError(char* lastSeenOption){
        cout<<_programb_<<" : "<<lastSeenOption<<" : Unknown option or incorrect syntax. Use -h for syntax summary."<<endl;
	exit(1);
}

void exitWithHelp()
{
	printHelp();
	//printCommandLineError();
	exit(1);
}

///Prints date with a leading # sign calling underlying OS. 
void systemDate(){
	system("printf \042#\042 ");
	system("date");
}

void printHelp(){
	cout<<"Usage:"<<endl;
	cout<<"       "<<_programb_<<" [-h|--help] (Use --help for more details)"<<endl;
	cout<<"       "<<_programb_<<" --version "<<endl;
	cout<<"       "<<_programb_<<" [OPTIONS] COMMAND ARGS "<<endl;
	cout<<endl;
	cout<<"        OPTIONS                                "<<endl;
	cout<<"                --debug                        "<<endl;
	cout<<"                --verbose                      "<<endl;
	cout<<"                -q , --quiet                   "<<endl;
	cout<<"                -z, --pid-normalization [s|p|r|l]"<<endl;
	cout<<"                -t , --format partition_format "<<endl;
	cout<<"                --tab tab_file                 "<<endl;
	cout<<"                --DIST_SUBSPROJECT             "<<endl;
	cout<<"                --beta beta_value              "<<endl;
	cout<<"                --mu mu_value                  "<<endl;
	cout<<endl;
	cout<<"        COMMANDS                               "<<endl;
	cout<<endl;
	cout<<"   For analyzing partitions                                          "<<endl;
        cout<<"       (-v|-e|-p|-i) partition1  partition2  [partition1_offset (=2) ] [partition2_offset (=partition1_offset) ] "<<endl;
	cout<<"       -c matrix-of-values partition1 [threshold (=-1.0)] [partition_offset (=2)]"<<endl;
	cout<<"       -d matrix-of-values partition1 [partition_offset (=2)]"<<endl;
	cout<<"       (-J|-R|-T) [-ext extensivity] [-ofs partition_offset (=2)] [-f partition_list | partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<"       (-V|-E|-P) [-ofs partition_offset (=2)] [-f partition_list | partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<"       --pstat-sym [-ofs partition_offset (=2)] [-f partition_list | partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<"       (--ipot|--cpot|--jpot|--v-measure-h) entropy [-ext extensivity] [-ofs partition_offset (=2)] [-f partition_list | partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<"       (-C|-H) [-cons] [-ofs partition_offset (=2)] [-f partition_list | [partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<"       (-A|-S|-I) [-ofs partition_offset (=2)] [-f partition_list | [partition1 [ partition2 [ ... ]]]"<<endl;
	cout<<endl;
	cout<<"   For editing partitions                                             "<<endl;
	cout<<"       --part-extract-elements elements_file [-tab mcl_tab_file] partition [partition1_offset (=2) ]" <<endl;
	cout<<"       --part-sort partition " <<endl;
	cout<<"       --part-sort-rename partition [prefix] " <<endl;
	cout<<"       --part-swap-names partition (requires use of --tab)" <<endl;
	cout<<endl;
	cout<<"   For converting between different partition formats                 "<<endl;
	cout<<"       --toMCL [-tab mcl_tab_file] partition [partition1_offset (=2) ]" <<endl;
	cout<<"       --toFREE partition [partition1_offset (=2) ]" <<endl;
	cout<<"       --MCLtoPART [-tab mcl_tab_file] partition [partition1_offset (=2) ]" <<endl;
	cout<<endl;
	cout<<"   For dealing with (fasta) sequence files                           "<<endl;
	cout<<"       --seq-noclone-sequences fasta_sequence_file [reference_sequence_file]"<<endl;
	cout<<endl;
	cout<<"   For analyzing Multiple Sequence Alignments                        "<<endl;
	cout<<"       --msa-seqid-stat [--positions positions_file] multiple_seq_alignment.fasta [multiple_seq_alignment.fasta2]"<<endl;
	cout<<"       --msa-seqid-avg [-thr threshold=50] multiple_seq_alignment.fasta [multiple_seq_alignment.fasta2]"<<endl;
	cout<<"       --msa-extract-positions positions_file multiple_seq_alignment.fasta"<<endl;
	cout<<"       --msa-extract-sequences sequences_file multiple_seq_alignment.fasta"<<endl;
	cout<<"       --msa-drop-sequences sequences_file multiple_seq_alignment.fasta"<<endl;
	cout<<"       --msa-extract-sequences-by-id  msa_file1 msa_file2 [minId maxId]"<<endl;
	cout<<"       --msa-drop-sequences-by-id sequences_file msa_file [minId maxId]"<<endl;
	cout<<"       --msa-extract-sequences-by-topid msa_file1 msa_file2 [count]"<<endl;
	cout<<"       --msa-drop-sequences-by-topid msa_file1 msa_file [count]"<<endl;
	cout<<"       --msa-map-partition partition multiple_seq_alignment.fasta "<<endl;
	cout<<"       --print-msa multiple_seq_alignment.fasta "<<endl;
	cout<<"       --msa-redundant [-nsam nsam] [-nseq nseq] [-seed seed] multiple_seq_alignment.fasta "<<endl;
	cout<<endl;
	cout<<"   For dealing with -interaction- matrices                           "<<endl;
	cout<<"       --edge-dist matrix-of-values [partition [partition_offset (=2)] ] "<<endl;
	cout<<"       -m matrix-of-values1 matrix-of-values2 "<<endl;
	cout<<"       -r matrix-of-values1 matrix-of-values2 partition [partition_offset (=2)]"<<endl;
	cout<<"       -l matrix-of-values1 matrix-of-values2 "<<endl;
	cout<<"       --print-matrix matrix-of-values "<<endl;
}
void printHelpLong(){
	cout<<"                                                                       "<<endl;
	cout<<_programb_<<" aims at being a general program for analyzing (sets of) partitions."<<endl;
	cout<<" Here a partition is defined as in set theory of mathematics (see      "<<endl;
	cout<<" http://en.wikipedia.org/wiki/Partition_of_a_set). It also allows to   "<<endl;
	cout<<" edit (rudimentarily), as well as generate, partitions.                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<" Whenever many input files are expected, one can either list them as   "<<endl;
	cout<<" command line arguments, or list them in a file and use option -f to   "<<endl;
	cout<<" specify that file.                                                    "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<" For calculating distances between partitions with different number of "<<endl;
	cout<<" elements, use option --DIST_SUBSPROJECT right before any *stat command."<<endl;
	cout<<" Works only with a *stat distance command, i.e., not  purity scores.   "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<" OPTIONS:                                                              "<<endl;
	cout<<"       -z , --pid-normalization norm                                   "<<endl;
	cout<<"          Determines the normalization used for calculating percent    "<<endl;
	cout<<"          sequence identities. The possible string values for norm are:"<<endl;
	cout<<"                 s , shorter-sequence                                  "<<endl;
	cout<<"                 p , aligned-positions                                 "<<endl;
	cout<<"                 r , aligned-residues                                  "<<endl;
	cout<<"                 l , average-length                                    "<<endl;
	cout<<"          Default normalization is the average sequence length, l.     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<" COMMANDS:                                                             "<<endl;
	cout<<"For analyzing partitions                                               "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -v , --vi-distance                                              "<<endl;
	cout<<"          Calculate VI distances between partition1 & partition2       "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -e , --edit-distance                                            "<<endl;
	cout<<"          Calculate the edit score distance between partition1 and     "<<endl;
	cout<<"               partition2                                              "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -p , -purity-scores                                             "<<endl;
	cout<<"          Calculates the purity scores of partition2 (the target)      "<<endl;
	cout<<"               againts the partition1 (the reference).                 "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -i , --intersection                                             "<<endl;
	cout<<"          Calculate the intersection of  partition1 & partition2       "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -c , --check-consistency-of-partition , --ccop [-tab tab_file]  "<<endl;
	cout<<"          Check cluster consistency according to the given matrix. If  "<<endl;
	cout<<"          partition and graph matrix label items differently, use the  "<<endl;
	cout<<"          option -tab to provide a tab file specifying the conversion. "<<endl;
	cout<<"          (See below for syntaxis of the matrix and tab file)          "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -d , --intra-inter-edge-dist                                    "<<endl;
	cout<<"          Calculate intra and inter cluster distribution of weights    "<<endl;
	cout<<"          according to the given matrix                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -J , --jstat [-ext extensivity] [-ref]                          "<<endl;
	cout<<"          Calculates Tarantola distance  for each pair of partitions.  "<<endl;
	cout<<"          For that it uses the Jeffrey's Qnorm based on Shannon Entropy."<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"          Default extensivity coefficient is "<<EXTENSIVITY_DEFAULT_RENYI<<"."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -R , --rstat [-ext extensivity] [-ref]                          "<<endl;
	cout<<"          Calculates Renyi distances for each pair of partitions.      "<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"          Default extensivity coefficient is "<<EXTENSIVITY_DEFAULT_RENYI<<"."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -T , --tstat [-ext extensivity] [-ref]                          "<<endl;
	cout<<"          Calculates Tsallis distances for each pair of partitions.    "<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"          Default extensivity coefficient is "<<EXTENSIVITY_DEFAULT_TSALLIS<<"."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -B , --bstat [-ref]                                             "<<endl;
	cout<<"          Calculates the Boltzmann distance for each pair of partitions."<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -V , --vstat [-ref]                                             "<<endl;
	cout<<"          Calculates the VI distance for each pair of partitions.      "<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -E , --estat [-ref]                                             "<<endl;
	cout<<"          Calculates the Edit Score distance for each pair of partitions"<<endl;
	cout<<"          With option -ref, the first partition is taken as a reference"<<endl;
	cout<<"          and it calculates the distances of all againts that one.     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -P , --pstat [-ref]                                             "<<endl;
	cout<<"          Calculates the purity scores (strict and lax) for each pair  "<<endl;
	cout<<"          of partitions. With option -ref, it calculates the purity    "<<endl;
	cout<<"          scores of all againts the first one, which is taken as a     "<<endl;
	cout<<"          reference.                                                   "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --pstat-sym , --pstat-symmetric                                 "<<endl;
	cout<<"          Calculates arithmetic averages of purity stric and purity lax"<<endl;
	cout<<"          scores for each pair of partitions.                          "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -n , --ipot  entropy [-ext extensivity]                         "<<endl;
	cout<<"          Calculates (information theoretic) potential (entropy) of each"<<endl;
	cout<<"          partition. The possible values for entropy are:              "<<endl;
	cout<<"          v | shannon                                                  "<<endl;
	cout<<"          e | cardinality                                              "<<endl;
	cout<<"          r | renyi                                                    "<<endl;
	cout<<"          t | tsallis                                                  "<<endl;
	cout<<"          Both, long and short option names are valid.                 "<<endl;
	cout<<"          Default extensivity coefficient is "<<EXTENSIVITY_DEFAULT_RENYI<<"."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --cpot, --conditional-potential entropy [-ext extensivity]      "<<endl;
	cout<<"          Calculates conditional entropy for each partition.           "<<endl;
	cout<<"          The possible values for entropy and extensivity are the same "<<endl;
	cout<<"          as for option --ipot.                                        "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --jpot, --joint-potential entropy [-ext extensivity]            "<<endl;
	cout<<"          Calculates joint entropy for each partition.                 "<<endl;
	cout<<"          The possible values for entropy and extensivity are the same "<<endl;
	cout<<"          as for option --ipot.                                        "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --v-measure-h , --v-measure-harmonic entropy [-ext extensivity] "<<endl;
	cout<<"          Calculates the Vmeasure between each pair of partitions. This"<<endl;
	cout<<"          measure is as that defined by Roseberg, A. and Hirschberg, J."<<endl;
	cout<<"          in http://acl.ldc.upenn.edu/D/D07/D07-1043.pdf. Use global   "<<endl;
	cout<<"          option --beta for specifying relative weight of homogeneity  "<<endl;
	cout<<"          versus completeness. Default is equal weight, i.e., beta=1.  "<<endl;
	cout<<"          and thus the average between both is strictly an harmonic one."<<endl;
	cout<<"          The possible values for entropy and extensivity are the same "<<endl;
	cout<<"          as for option --ipot.                                        "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --v-measure-a , --v-measure-arithmetic entropy [-ext extensivity]"<<endl;
	cout<<"          Analogous to --v-measure-h but using arithmetic mean between "<<endl;
	cout<<"          homogeneity and completeness.                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --v-measure-g , --v-measure-geometric entropy [-ext extensivity]"<<endl;
	cout<<"          Analogous to --v-measure-h but using geometric mean between "<<endl;
	cout<<"          homogeneity and completeness.                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -C , --cluster-stat [-ofs ofs] [-norm gaug] [-cons|-consensus]  "<<endl;
	cout<<"          For each item, determines the most frequent cluster where    "<<endl;
	cout<<"          it appears among all the clusters of all the given partitions."<<endl;
	cout<<"          It also prints its size and observed frequency (both, raw    "<<endl;
	cout<<"          count and %).                                                "<<endl;
	cout<<"          Option -ofs,see below, allows to specify a partition offset. "<<endl;
	cout<<"          Option -norm gaug gauges the normalization used for          "<<endl;
	cout<<"                 determining the %frequencies. By default these are    "<<endl;
	cout<<"                 calculated by counting how many times the mode cluster"<<endl;
	cout<<"                 is found at each of the different partitions and then "<<endl;
	cout<<"                 dividing by the number of partitions N. With this     "<<endl;
	cout<<"                 option, that count gets divided by N+gaug, where gaug "<<endl;
	cout<<"                 can be negative or positive.                          "<<endl;
	cout<<"          Option -cons or -consensus will print the consensus partition"<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -A , --adjacency-stat , {--adjstat}                             "<<endl;
	cout<<"          Determines the average adjacency matrix from the provided    "<<endl;
	cout<<"          partitions. The adjacency matrix of a partition is the graph "<<endl;
	cout<<"		 where edges (0 or 1 ) represent two elements belonging to the"<<endl;
	cout<<"		 same subfamily. The average adjacency matrix has edges with  "<<endl;
	cout<<"		 continous values [0,1]. The output consists in a matrix of   "<<endl;
	cout<<"		 values and a gray-scale image of it in PGM format.           "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -S , --split-merge-analysis , {--splitstat}                     "<<endl;
	cout<<"         (Split-Merge plot)                                            "<<endl;
	cout<<"          Determines the overlap of each cluster to those of the       "<<endl;
	cout<<"          reference partition (the first). Possible values are for     "<<endl;
	cout<<"          the overlap are:                                             "<<endl;
	cout<<"          -over fraction elements in common relative to the target cluster."<<endl;
	cout<<"          -cos  cosine normalized similarity                           "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"          It outputs:                                                  "<<endl;
	cout<<"          -Confusion matrix (in % of the target clusters) taking       "<<endl;
	cout<<"            the first partition as reference and the second as target. "<<endl;
	cout<<"          -number of overlaps for each target cluster                  "<<endl;
	cout<<"          -Split-Merge image showing the CT. In addition it show two   "<<endl;
	cout<<"		   reference color bars: a bottom color bar representing the  "<<endl;
	cout<<"            perfect split transformations (black), the merge-only      "<<endl;
	cout<<"            (white) ones and those cases in between (different grey    "<<endl;
	cout<<"            levels); a right-most column shows whether these are perfect"<<endl;
	cout<<"            matches (black) or not (white).                            "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -I , --info , {--isPart , --isaPart , --is-partition} [-ofs ofs]"<<endl;
	cout<<"           For each partition checks whether it is a sound partition   "<<endl;
	cout<<"           or not, i.e., whether all of its clusters are pair-wise     "<<endl;
	cout<<"           disjoint. With option -q, only error message will be printed"<<endl;
	cout<<"           in case partition is not sound, otherwise it'll keep silent."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -H, --hasse-diagram                                             "<<endl;
	cout<<"           prints the local Hasse Diagram (graph) spanned by the       "<<endl;
	cout<<"           given partitions.                                           "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"For editing partitions                                                 "<<endl;
	cout<<"       --part-extract-elements {--extract-elements} elements_file      " <<endl;
	cout<<"           elements_file lists the names of the elements to cull from  "<<endl;
	cout<<"           the given partition                                         "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --part-sort                                                     "<<endl;
	cout<<"           Sorts the clusters by size, the larger on top. Ties are     "<<endl;
	cout<<"           sorted alphabetically by their first item. Within each      "<<endl;
	cout<<"           cluster, items are sorted alphabetically.                   "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --part-sort-rename partition [prefix]                           "<<endl;
	cout<<"           As --part-sort, but also rename each clusters consecutively "<<endl;
	cout<<"           as C1, C2,etc. If a prefix string is supplied use that      "<<endl;
	cout<<"           instead of C.                                               "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --part-swap-names partition                                     "<<endl;
	cout<<"       --part-swap-labels partition                                    "<<endl;
	cout<<"           Swaps elements' names present in partition by their new     "<<endl;
	cout<<"           names as found in the provided tab file. An element's name in "<<endl;
	cout<<"           the partition will be changed iif there is a translation for"<<endl;
	cout<<"           it found in the tab file; otherwise it will be left as it is."<<endl;
	cout<<"           Thus, it is not mandatory to provide a translation for all  "<<endl;
	cout<<"           elements. Requires the use of --tab to specify a tab file   "<<endl;
	cout<<"           providing the mapping between new and old names. See general"<<endl;
	cout<<"           options.                                                    "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"For converting between different partition formats                     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --toMCL [-tab mcl_tab_file] converts partition from PART format  "<<endl;
	cout<<"           to MCL's format. If additional tab file is provided, output "<<endl;
	cout<<"           will contain the specific label index given in the tab file."<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --toFREE converts partition from PART format to FREE format.     "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --MCLtoPART [-tab mcl_tabl_file]                                 "<<endl;
	cout<<"           converts partition from MCL format to PART format.          "<<endl;
	cout<<"           If additional tab file is provided, output will contain     "<<endl;
	cout<<"           items' labels, instead of simply their MCL index number.    "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"For dealing with (fasta) sequence files                                "<<endl;
	cout<<"       --drop-clone-sequences                                          "<<endl;
	cout<<"       --msa-noclone-sequences                                         "<<endl;
	cout<<"       --seq-noclone-sequences sequence_file (fasta)                   "<<endl;
	cout<<"          Given a (fasta) sequence file or a fasta MSA file, remove all"<<endl;
	cout<<"          duplicate sequences. Here duplicate means literally that,    "<<endl;
	cout<<"          namely, exactly the same string of characters. Therefore, it "<<endl;
	cout<<"          is not the same as having a pid=100%, but more stringent.    "<<endl;
	cout<<"          If a second sequence file is provided, drop also sequences   "<<endl;
	cout<<"          that are clones of any sequence in the second file.          "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"For analyzing Multiple Sequence Alignments                             "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-seqid-stat                                                "<<endl;
	cout<<"       --msa-seqid-stat [--positions file]                             "<<endl;
	cout<<"          Given a multiple sequence alignment in fasta format, it      "<<endl;
	cout<<"          prints all pair-wise sequence identities. By default, it     "<<endl;
	cout<<"          calculates identities over the full sequence length. The     "<<endl;
	cout<<"          second version allows to specify the (reduced) set of positions"<<endl;
	cout<<"          we want to consider in comparing sequences. These should be  "<<endl;
	cout<<"          specified in a file, each separated by space,tabs, new lines,"<<endl;
	cout<<"          etc. The positions are understood as columns of the MSA.     "<<endl;
	cout<<"          If two MSA are provided, it prints the sequence Id of the    "<<endl;
	cout<<"          first set against the second.                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-seqid-avg [-thr threshold ]                               "<<endl;
	cout<<"          Similar as option --msa-seqid-avg, but prints for each sequence"<<endl;
	cout<<"          a statistics of its pair-wise sequence identity to all other "<<endl;
	cout<<"          sequences. This consists of average Seq.Id, standard         "<<endl;
	cout<<"          deviation, variance, minimum Seq.Id, maximum Seq.Id, number  "<<endl;
	cout<<"          of pairs with Seq.Id > threshold, fraction of pairs with Seq."<<endl;
	cout<<"          Id. > threshold and total number of pairs.                   "<<endl;
	cout<<"          Option -thr allows to provide a specific threshold to use.   "<<endl;
	cout<<"                      default value is 50%. Values are floating numbers"<<endl;
	cout<<"                      within [0,100].                                  "<<endl;
	cout<<"          If two MSA are provided, it prints the sequence Id of the    "<<endl;
	cout<<"          first set against the second.                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-extract-positions positions_file msa_file                 "<<endl;
	cout<<"          From the given MSA, extract only columns specified in file   "<<endl;
	cout<<"          positions_file.                                              "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-extract-sequences sequences_file msa_file                 "<<endl;
	cout<<"       --msa-drop-sequences sequences_file msa_file                    "<<endl;
	cout<<"          From the given MSA, extract only sequences specified in file "<<endl;
	cout<<"          sequences_file. This file contains a list of sequences names "<<endl;
	cout<<"          The second form drops those sequences instead.               "<<endl;
	cout<<"          If a positions file is given, sequence Id's are calculated   "<<endl;
	cout<<"          considering only those columns of the MSA.                   "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-extract-sequences-by-id msa_file1 msa_file2 [minId maxId] "<<endl;
	cout<<"       --msa-drop-sequences-by-id sequences_file msa_file [minId maxId]"<<endl;
	cout<<"          From MSA msa_file1, extract sequences with an ID above minId "<<endl;
	cout<<"          and at most maxId against any sequence of MSA msa_file2.     "<<endl;
	cout<<"          The second form drops those sequences instead. Default values"<<endl;
	cout<<"          values are minId=30 and maxId=100, i.e., homologous sequences."<<endl;
	cout<<"          If a positions file is given, sequence Id's are calculated   "<<endl;
	cout<<"          considering only those columns of the MSA. In this case minId"<<endl;
	cout<<"          and maxId are mandatory and must come before positions_file. "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-extract-sequences-by-topid msa_file1 msa_file2 [count]    "<<endl;
	cout<<"       --msa-drop-sequences-by-topid sequences_file msa_file [count]   "<<endl;
	cout<<"          From MSA msa_file1, extract at most count most similar sequences "<<endl;
	cout<<"          (seq.ID) to any sequence of MSA msa_file2.                   "<<endl;
	cout<<"          The second form drops those sequences instead.               "<<endl;
	cout<<"          If a positions file is given, sequence Id's are calculated   "<<endl;
	cout<<"          considering only those columns of the MSA. In this case count"<<endl;
	cout<<"          is mandatory and must come before positions_file.            "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-redundant [-nsam nsam] [-nseq nseq] [-seed seed]          "<<endl;
	cout<<"          Duplicates sequences chosen at random in the given multiple  "<<endl;
	cout<<"          sequence alignment. Wtihout options, only one is chosen.     "<<endl;
	cout<<"          Option -nsam  Generate nsam samples of MSAs with nseq dupli- "<<endl;
	cout<<"                        cated sequences. Each sample is written is its "<<endl;
	cout<<"                        own directory.                                 "<<endl;
	cout<<"                 -nseq  Specify the number of sequences to duplicate.  "<<endl;
	cout<<"                 -seed  Specify the seed of the random number generator"<<endl;
	cout<<"          All options are expected to be integer values. The value of  "<<endl;
	cout<<"          the seed is written within .seed_used allowing for repeated  "<<endl;
	cout<<"          experiments.                                                 "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --msa-map-partition                                             "<<endl;
	cout<<"          Given a Partition and the original MSA, output the MSA       "<<endl;
	cout<<"          with the cluster annotation of the SDPpred server            "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --print-msa                                                     "<<endl;
	cout<<"          Prints the given multiple sequence alignment. Useful for     "<<endl;
	cout<<"          debugging.                                                   "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"For dealing with -interaction- matrices                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --edge-dist                                                     "<<endl;
	cout<<"          For each node, prints the distribution of edge weights.      "<<endl;
	cout<<"          Information printed is: Node, average edge weight, standard  "<<endl;
	cout<<"          deviation, standard error, skewness, minimum edge value,     "<<endl;
	cout<<"          max edge value and sample size (number of edges).            "<<endl;
	cout<<"          If a partition is provided, it also prints the cluster size  "<<endl;
	cout<<"          and cluster name each node belongs to.                       "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -m , --merge-graphs                                             "<<endl;
	cout<<"          Merge two graph matrices into one that contains both values  "<<endl;
	cout<<"          for each pair of items, i.e., the resulting graph looks like "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"                        stringA stringB  float1 float2                 "<<endl;
	cout<<"                           ...    ...      ...    ...                  "<<endl;
	cout<<"          where float1, float2 are the matrix values of matrix1 and    "<<endl;
	cout<<"          matrix2, respectively. Both matrices are expected to contain "<<endl;
	cout<<"          the same set of pair of items, i.e., the same set of edges.  "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -r , --merge-graphs-color                                       "<<endl;
	cout<<"          as option -m, but in addition includes the name of the   "<<endl;
	cout<<"          cluster each pair of values belong to. If they belong to     "<<endl;
	cout<<"          different clusters the label is \042x\042. The label is NAN  "<<endl;
	cout<<"          if any of the item does not belong to any of the clusters    "<<endl;
	cout<<"          defined in the given partition. The format of the output is  "<<endl;
	cout<<"                 float1A flaot2 clustername_AB stringA stringB         "<<endl;
	cout<<"                   ...     ...      ...          ...    ...            "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -l , --cull-edges                                               "<<endl;
	cout<<"          Culls from matrix of values edges specified in second file.  "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --print-matrix                                                  "<<endl;
	cout<<"          Print the given interaction matrix. For debugging.       "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"General options                                                        "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       --verbose                                                       "<<endl;
	cout<<"          For debugging.                                               "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -q , --quiet                                                    "<<endl;
	cout<<"          quiet mode. Do not print out comment lines (that start with  "<<endl;
	cout<<"          `#').                                                        "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       -t , --format {--fmt} [pfmt=input_partition_file_format]        "<<endl;
	cout<<"          Possible format values are: PART,MCL and FREE. See below.    "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"File formats:                                                          "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       matrix-of-values:                                               "<<endl;
	cout<<"                        stringA stringB  float                         "<<endl;
	cout<<"                        stringA stringC  float                         "<<endl;
	cout<<"                          ...     ...     ...                          "<<endl;
	cout<<"                        stringZ stringV  float                         "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       tab file:                                                       "<<endl;
	cout<<"                        integer1  string1                              "<<endl;
	cout<<"                        integer2  string2                              "<<endl;
	cout<<"                          ...       ...                                "<<endl;
	cout<<"                                                                       "<<endl;
	cout<<"       partition:                                                      "<<endl;
	cout<<"         PART: (default, i.e., partition_offset=2)                     "<<endl;
	cout<<"               sizeA clusterA_name  item_1 item_2 ... item_sizeA       "<<endl;
	cout<<"               sizeB clusterB_name  item_1 item_2 ... item_sizeB       "<<endl;
	cout<<"                ...      ...       ...    ...         ...              "<<endl;
	cout<<"            or (partition_offset=1):                                   "<<endl;
	cout<<"               sizeA  item_1 item_2 ... item_sizeA                     "<<endl;
	cout<<"                ...    ...    ...         ...                          "<<endl;
	cout<<"         FREE: (not yet implemented) (partition_offset=0)              "<<endl;
	cout<<"               item_1 item_2 ... item_sizeA                            "<<endl;
	cout<<"                ...    ...         ...                                 "<<endl;
	cout<<"         MCL : MCL's own matrix format for partitions. See MCL manual. "<<endl;
}

#endif //END _HELPFUNCTIONS
