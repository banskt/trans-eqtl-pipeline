## Format output as: ensembl_id	chrom	start	end	name
BEGIN { FS=OFS="\t" }
{
  sub("chr", "", $1)
  split($4, meta, ";")
  
  if ($1 != "X" && $1 != "Y" && $1 != "M") {
    for (i in meta) {
      ## remove leading whitespace
      gsub (/^ */, "", meta[i])
      ## find Ensembl ID and put it in a variable
      if (meta[i] ~ /^piRNA_code/) {
        split(meta[i], arr01, " ") # split with whitespace
        gsub(/\"/, "", arr01[2]) # remove quotes from second element
        ensid_trim = arr01[2]
      }
    }

    print ensid_trim, $1, $2, $3;
  }
}
