#
# --- Test badnex.R ---
#

if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexusfiles")
} else {
    pth <- system.file(package="phylobase", "nexusfiles")
}

badFile <- file.path(pth, "badnex.nex")

test.checkTree <- function() {
    checkException(readNexus(file=badFile))
}

