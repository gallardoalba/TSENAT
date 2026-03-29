# Replace Nullable<int> n_tx_fixed = R_NilValue with int n_tx_fixed = -1
s/Nullable<int> n_tx_fixed = R_NilValue/int n_tx_fixed = -1/g
s/Nullable<int> n_transcripts_fixed = R_NilValue/int n_transcripts_fixed = -1/g

# Replace isNotNull checks with > 0 checks
s/n_tx_fixed\.isNotNull()/n_tx_fixed > 0/g
s/n_transcripts_fixed\.isNotNull()/n_transcripts_fixed > 0/g

# Replace as<int> conversions to just use the value directly
s/n_tx_use = as<int>(n_tx_fixed);/n_tx_use = n_tx_fixed;/g
s/n_tx_use = as<int>(n_transcripts_fixed);/n_tx_use = n_transcripts_fixed;/g
