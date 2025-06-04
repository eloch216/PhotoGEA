# The debug output includes timestamps, which will never be identical between
# runs. Some of the quotes are altered when saving/reading stored outputs, which
# can also cause problems.
process_saved_debug <- function(debug_output) {
    # Remove timestamps
    debug_output <- gsub(
        '^Time: [0123456789-]+ [0123456789:.]+ +',
        '',
        debug_output
    )

    # Standardize quotes
    gsub(
        '[‘’]',
        "'",
        debug_output
    )
}
