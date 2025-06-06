#!/bin/bash

# Usage: download from google drive
#   1. Processed input files for pgenmi and fwpgenmi runs 
#   2. Raw histone marks peaks, ATACseq, DE gene p-values for generating input files


# before uploading to google drive 
# cd {PROJPATH}/input
# tar -czf evidences.tar.gz TFBS_*
# https://drive.google.com/file/d/1saLg9f93Bne0dNHAUWXjZwuUolF9pcoE/view?usp=share_link
# 
# do the same with raw data
# cd {PROJPATH}/preprocess/data
# tar -czf rawdata.tar.gz ATAcSeq  deseq2  K27ac  K27me3  K4me1  K4me3
# https://drive.google.com/file/d/1hnqe0hfHRoHF_01uU8jhP0gFLFnCjIy8/view?usp=share_link
#
# Uploading to google drive:
# https://stackoverflow.com/questions/58589734/pydrive-trying-to-upload-files-to-google-drive-from-a-remote-server
# https://pythonhosted.org/PyDrive/oauth.html

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
PROJPATH=$SCRIPT_PATH/..
source $PROJPATH/activate venv

# -------------------------------------------------------------------------- # 
#                   Get the preprocessed evidence files for                  #
#  TFBS_DiffMark, TFBS_DiffACC, TFBS_PresMark, TFBS_PresACC, and TFBS_only   #
# -------------------------------------------------------------------------- #
echo "Downloading the preprocess evidence file from google drive ..."
mkdir -p $PROJPATH/input
cd $PROJPATH/input
EVIDID="1saLg9f93Bne0dNHAUWXjZwuUolF9pcoE"
gdown $EVIDID
# sanity check
#if ! md5sum -c input.md5; then echo "corrupted files"; fi
tar -xzf evidences.tar.gz
rm evidences.tar.gz

# ------------------------------------------------------------------------------- # 
#    Get the multi-omics data containing histone mark peaks, ATACseq, DE genes    #
# ------------------------------------------------------------------------------- #
echo "Downloading the multi-omics data for preprocessing from google drive ..."
mkdir -p $PROJPATH/preprocess/data
cd $PROJPATH/preprocess/data
RAWID="1hnqe0hfHRoHF_01uU8jhP0gFLFnCjIy8"
gdown $RAWID
# sanity check
#if ! md5sum -c input.md5; then echo "corrupted files"; fi
tar -xzf rawdata.tar.gz
rm rawdata.tar.gz

# -------------------------------------------------- # 
#     Output files of some preprocessing stages      #
# -------------------------------------------------- #
echo "Downloading output files for some preprocessing stages ..."
cd $PROJPATH/preprocess
RAWID="1n6GTxuXaMlqP7icDha_JF_9dYL3reQ6e"
gdown $RAWID
tar -xzf preprocess_stages_outfiles.tar.gz
rm preprocess_stages_outfiles.tar.gz