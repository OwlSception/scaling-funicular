setup_logger <- function(logfile = "app.log", level = "INFO") {
    if (!requireNamespace("logger", quietly = TRUE)) {
        install.packages("logger")
    }
    library(logger)
    
    # Set log layout
    log_layout(layout_glue_generator(format = '[{level}] {time} {msg}'))
    
    # Set log threshold
    log_threshold(level)
    
    # Set log appender to file
    log_appender(appender_tee(logfile))
    
    invisible(TRUE)
}

