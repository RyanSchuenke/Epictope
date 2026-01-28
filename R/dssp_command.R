#' @title dssp_command
#'
#' @description Runs the DSSP command with mkdssp to process a mmCIF file and generates a DSSP file.
#'
#' @param cif_file A character vector specifying the path to the input mmCIF file.
#'
#' @return The function returns the name of the output DSSP file.
#'
#' @examples
#' # Generate the DSSP file from the input mmCIF file "protein.cif"
#' #dssp_command("protein.cif")
#'
#' @export

dssp_command <- function(cif_file) {
  # Check that .x is a character vector
  stopifnot(is.character(cif_file))

  # Define output filename with config model_folder and replacing the extension with ".dssp"
  output_file <- file.path(model_folder, basename(gsub("\\.cif.gz|\\.cif", ".dssp", cif_file)))

  # Run the DSSP command with mkdssp, return error message if an error occurs.
  system(command = paste0("mkdssp ", cif_file, " ", output_file), intern = TRUE)

  # Return the name of the output file.
  return(output_file)
}
