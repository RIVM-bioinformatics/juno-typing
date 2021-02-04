#!/usr/bin/env bats

## General Juno tests

@test "Printing help" {
  printed_help=`bash juno-typing.sh -h`
  [[ "${printed_help}" =~ "-i, --input [DIR]                 This is the folder containing your input" ]]
}

@test "Make sample sheet from input directory that does not contain either fasta or fastq files" {
  python bin/generate_sample_sheet.py tests/ > tests/test_sample_sheet.yaml
  sample_sheet_errors=`diff --suppress-common-lines tests/test_sample_sheet.yaml tests/example_output/wrong_sample_sheet.yaml`
  [[ -z $sample_sheet_errors ]]
}

@test "Make sample sheet from fastq or mixed input" {
  python bin/generate_sample_sheet.py tests/example_fastq_input/ > tests/test_sample_sheet.yaml
  sample_sheet_errors=`diff --suppress-common-lines tests/test_sample_sheet.yaml tests/example_output/fastq_sample_sheet.yaml`
  [[ -z $sample_sheet_errors ]]
}

@test "Make sample sheet from fasta input" {
  python bin/generate_sample_sheet.py tests/example_fasta_input/ > tests/test_sample_sheet.yaml
  sample_sheet_errors=`diff --suppress-common-lines tests/test_sample_sheet.yaml tests/example_output/fasta_sample_sheet.yaml`
  [[ -z $sample_sheet_errors ]]
}

## Specific for MLST7

@test "Making sample sheet fastq + metadata (species name)" {
  python bin/generate_sample_sheet.py tests/example_fastq_input --metadata tests/files/example_metadata.csv > tests/test_sample_sheet.yaml
  sample_sheet_errors=`diff --suppress-common-lines tests/test_sample_sheet.yaml tests/example_output/metadata_sample_sheet.yaml`
  [[ $(cat tests/test_sample_sheet.yaml) =~ "senterica" ]]
  [[ -z $sample_sheet_errors ]]
}

@test "Test full pipeline (dry run)" {
  bash juno-typing.sh -i tests/example_fastq_input/ --species Salmonella enterica -n
  [[ "$status" -eq 0 ]]
}

@test "Test error occurs when neither species nor metadata file are provided" {
  run juno-typing.sh -i tests/example_fastq_input/ -n
  [[ ! "$status" -eq 0 ]]
}

@test "Check full pipeline if running locally (and test samples present)" {
  bash juno-typing.sh -i tests/example_fastq_input/ --metadata tests/files/example_metadata.csv
  [[ "$status" -eq 0 ]]
}
