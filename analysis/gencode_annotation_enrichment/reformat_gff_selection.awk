## Format output as: ensembl_id	chrom	start	end	name
BEGIN { FS=OFS="\t" }
{
  sub("chr", "", $1)
  split($4, meta, ";")
  
  if ($1 != "X" && $1 != "Y" && $1 != "M") {
    for (i in meta) {
      ## remove leading whitespace
      gsub (/^ */, "", meta[i])
      ## find gene_name and put it in a variable
      if (meta[i] ~ /^Name/) {
        split(meta[i], arr01, "=") # split with '='
        gsub(/\"/, "", arr01[2]) # remove quotes from second element
        name = arr01[2]
      }
      ## find Ensembl ID and put it in a variable
      if (meta[i] ~ /^ID/) {
        split(meta[i], arr01, "=") # split with '='
        gsub(/\"/, "", arr01[2]) # remove quotes from second element
        #split(arr01[2], arr02, ".") 
        #ensid_full = arr01[2] 
        #ensid_trim = arr02[1]
        ensid_trim = arr01[2]
      }
    }

    print ensid_trim, $1, $2, $3, name;
  }
}
