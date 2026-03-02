#!/bin/bash

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc',).csv" --denovo-output helaqc_results/helaqc_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf --db-output helaqc_results/helaqc_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.10',).csv" --denovo-output helaqc_results/helaqc_decoy_0.10_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.10.mgf --db-output helaqc_results/helaqc_decoy_0.10_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.20',).csv" --denovo-output helaqc_results/helaqc_decoy_0.20_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.20.mgf --db-output helaqc_results/helaqc_decoy_0.20_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.30',).csv" --denovo-output helaqc_results/helaqc_decoy_0.30_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.30.mgf --db-output helaqc_results/helaqc_decoy_0.30_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.40',).csv" --denovo-output helaqc_results/helaqc_decoy_0.40_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.40.mgf --db-output helaqc_results/helaqc_decoy_0.40_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.50',).csv" --denovo-output helaqc_results/helaqc_decoy_0.50_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.50.mgf --db-output helaqc_results/helaqc_decoy_0.50_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.60',).csv" --denovo-output helaqc_results/helaqc_decoy_0.60_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.60.mgf --db-output helaqc_results/helaqc_decoy_0.60_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.70',).csv" --denovo-output helaqc_results/helaqc_decoy_0.70_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.70.mgf --db-output helaqc_results/helaqc_decoy_0.70_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.80',).csv" --denovo-output helaqc_results/helaqc_decoy_0.80_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.80.mgf --db-output helaqc_results/helaqc_decoy_0.80_database_results.csv

novoboard preprocess --denovo-file "/home/j-daniel/Documents/winnow/novoboard/preds_('helaqc.mgf.decoy_0.90',).csv" --denovo-output helaqc_results/helaqc_decoy_0.90_de_novo_results.csv --db-mgf-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf.decoy_0.90.mgf --db-output helaqc_results/helaqc_decoy_0.90_database_results.csv
