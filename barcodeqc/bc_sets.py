"""Barcode set dictionary used by merChecker and tests.

Each entry contains a `title`, a `merPath_template` which should be
expanded with `.format(merDir=merDir)` after `merDir` is defined in
the main script, and lists for `platesL1`, `platesL2` and `plateList`.
"""

bc_sets = {
    "bc50": {
        "title": "P50 Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList50.tsv",
        "platesL1": ["BCA"],
        "platesL2": ["BCB"],
        "plateList": [["BCA"], ["BCB"], ["BCB"]],
    },
    "bc96": {
        "title": "P96 Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList96.tsv",
        "platesL1": ["BCA_R_96", "BCA_L_96"],
        "platesL2": ["BCB_R_96", "BCB_L_96"],
        "plateList": [["BCA_R_96", "BCA_L_96"],
                      ["BCB_R_96", "BCB_L_96"],
                      ["BCB_R_96", "BCB_L_96"]],
    },
    "fg96": {
        "title": "FG96 Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merListfg96.tsv",
        "platesL1": ["BCAA_FG_96"],
        "platesL2": ["BCBB_FG_96"],
        "plateList": [["BCAA_FG_96"],
                      ["BCBB_FG_96"],
                      ["BCBB_FG_96"]],
    },
    "bc210": {
        "title": "P210 Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList210.tsv",
        "platesL1": ["BCA01", "BCA02", "BCA03", "BCA04"],
        "platesL2": ["BCB01", "BCB02", "BCB03", "BCB04"],
        "plateList": [["BCA01", "BCA02", "BCA03", "BCA04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"]],
    },
    "fg210p30": {
        "title": "P210p30 Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList210p30Draft.tsv",
        "platesL1": ["BCA01", "BCA02", "BCA03", "BCA04", "BCA05"],
        "platesL2": ["BCB01", "BCB02", "BCB03", "BCB04", "BCB05"],
        "plateList": [["BCA01", "BCA02", "BCA03", "BCA04", "BCA05"],
                      ["BCB01", "BCB02", "BCB03", "BCB04", "BCB05"],
                      ["BCB01", "BCB02", "BCB03", "BCB04", "BCB05"]],
    },
    "bc220_05-OCT": {
        "title": "220bc 05-OCT Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList220_05-OCT.tsv",
        "platesL1": ["BCA201", "BCA202", "BCA203"],
        "platesL2": ["BCB201", "BCB202", "BCB203"],
        "plateList": [["BCA201", "BCA202", "BCA203"],
                      ["BCB201", "BCB202", "BCB203"],
                      ["BCB201", "BCB202", "BCB203"]],
    },
    "bc220": {
        "title": "220bc 25-APR Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList220_25-APR.tsv",
        "platesL1": ["BCA201", "BCA202", "BCA203"],
        "platesL2": ["BCB201", "BCB202", "BCB203"],
        "plateList": [["BCA201", "BCA202", "BCA203"],
                      ["BCB201", "BCB202", "BCB203"],
                      ["BCB201", "BCB202", "BCB203"]],
    },
    "bc220_18-SEP": {
        "title": "220bc 18-SEP Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList220_18-SEP.tsv",
        "platesL1": ["BCA01", "BCA02", "BCA03", "BCA04"],
        "platesL2": ["BCB01", "BCB02", "BCB03", "BCB04"],
        "plateList": [["BCA01", "BCA02", "BCA03", "BCA04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"]],
    },
    "bc220_20-MAY": {
        "title": "220bc 20-MAY Barcode Set [bcSet: {bcSet}]",
        "merPath_template": "{merDir}/merList220_20-MAY.3.tsv",
        "platesL1": ["BCA01", "BCA02", "BCA03", "BCA04"],
        "platesL2": ["BCB01", "BCB02", "BCB03", "BCB04"],
        "plateList": [["BCA01", "BCA02", "BCA03", "BCA04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"],
                      ["BCB01", "BCB02", "BCB03", "BCB04"]],
    },
}

__all__ = ["bc_sets"]
