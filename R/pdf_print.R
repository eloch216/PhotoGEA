pdf_print <- function(
    plot_obj,
    width = 7,
    height = 7,
    save_to_pdf = FALSE,
    file = NULL,
    new_window = TRUE,
    ...
)
{
    if (save_to_pdf) {
        if (is.null(file)) {
            grDevices::pdf(width = width, height = height, ...)
        } else {
            grDevices::pdf(file = file, width = width, height = height, ...)
        }
        print(plot_obj)
        grDevices::dev.off()
    } else {
        if (new_window) {
            grDevices::dev.new(width = width, height = height)
        }
        print(plot_obj)
    }
}
