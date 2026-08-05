## ================================================================
## extract_helping_task_session.R
##
## Parses ONE raw Helping Task session file straight from the
## Pavlovia/jsPsych platform export -- the giant CSV-with-embedded-
## escaped-JSON file that looks like:
##   "trial_type","trial_index","time_elapsed","internal_node_id","value"
##   ...
##       ""trial_history_columns"": ""overall_idx,block_num,...""
##       ""trial_history"": [
##           ""1;4;1;1782392856739;...;4"",
##           ""2;4;2;...."",
##           ...
##       ],
##
## -- into the same clean column set used throughout this project:
## PID, Task Design, Trial, Block, Helper Policy, Effort Intended,
## Effort Actual, 3 reaction times, Satisfaction Rating, Trust Rating,
## Points For Helper, Points For Help-Seeker.
##
## ================================================================
## HOW TO USE
## 1. Set Working Directory to Data File Location
## 2. Scroll to Bottom and Update values at the bottom in "df <- ..."
## 3. Update File Name in "write.csv(...)"
## 4. Run Code
## ================================================================

suppressMessages(library(dplyr))

effort_numeric_map <- c(
  "Minimum" = 1,
  "1/3"     = 2,
  "2/3"     = 3,
  "Maximum" = 4
)

policy_map <- c("C" = "CARING", "U" = "UNCARING")

parse_helping_task_session <- function(filepath, pid, task_design,
                                          effort_map = effort_numeric_map,
                                          policy_lookup = policy_map,
                                          renumber_blocks = FALSE) {

  lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")

  ## ---- 1. get the column names from the "trial_history_columns" line ----
  header_line <- grep("trial_history_columns", lines, value = TRUE)[1]
  if (is.na(header_line)) {
    stop("Could not find a 'trial_history_columns' line in ", filepath,
         ". Is this a raw Helping Task session export?")
  }
  header_str <- sub('.*trial_history_columns"":\\s*""([^"]*)"".*', '\\1', header_line)
  col_names  <- strsplit(header_str, ",")[[1]]
  n_cols <- length(col_names)
  cat("Found", n_cols, "columns:", paste(col_names, collapse = ", "), "\n")

  ## ---- 2. pull out every trial row --------------------------------------
  ## Every trial line has the form   ""<idx>;<field>;<field>;...""[,]
  ## We first grab every line that LOOKS like it starts that way, then keep
  ## only the ones with exactly the right number of semicolons (n_cols-1) --
  ## this is what correctly excludes the platform's "events" log lines,
  ## which have the same "<number>;<text>" shape but only 1 semicolon.
  candidate_lines <- grep('^\\s*""[0-9]+;', lines, value = TRUE)
  semi_counts <- lengths(regmatches(candidate_lines, gregexpr(";", candidate_lines)))
  trial_lines <- candidate_lines[semi_counts == (n_cols - 1)]
  if (length(trial_lines) == 0) {
    stop("Found the header (", n_cols, " columns) but no trial rows with a ",
         "matching number of fields in ", filepath, ". Check the file wasn't ",
         "altered/re-saved by Excel (which can change the quoting).")
  }

  split_one <- function(x) {
    inner <- sub('^\\s*""(.*)""\\s*,?\\s*$', '\\1', x)
    strsplit(inner, ";")[[1]]
  }
  parsed <- lapply(trial_lines, split_one)

  ok_length <- sapply(parsed, length) == n_cols
  if (any(!ok_length)) {
    warning(sum(!ok_length), " trial row(s) did not have the expected ",
            n_cols, " fields and were dropped. Inspect trial_lines[!ok_length] ",
            "if this seems wrong.")
  }
  parsed <- parsed[ok_length]

  d <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
  names(d) <- col_names
  cat("Parsed", nrow(d), "trials.\n")

  ## ---- 3. numeric conversions --------------------------------------------
  numeric_cols <- c("overall_idx", "block_num", "block_idx",
                      "effort_selection_rxn_time", "effort_points",
                      "partner_points", "rate_self_rxn_time", "rate_self_value",
                      "rate_partner_rxn_time", "rate_partner_value")
  for (cc in numeric_cols) {
    if (cc %in% names(d)) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  }

  ## ---- 4. recode effort labels & helper policy ---------------------------
  recode_effort <- function(x) {
    if (is.null(effort_map)) return(x)
    out <- unname(effort_map[x])
    unmapped <- unique(x[is.na(out) & !is.na(x)])
    if (length(unmapped) > 0) {
      warning("Effort label(s) not found in effort_numeric_map: ",
              paste(unmapped, collapse = ", "),
              " -- these became NA. Add them to effort_numeric_map at the ",
              "top of this script and re-run.")
    }
    out
  }
  recode_policy <- function(x) {
    out <- unname(policy_lookup[x])
    out[is.na(out)] <- x[is.na(out)]  # keep raw code if not in lookup
    out
  }

  ## ---- 5. assemble exactly the requested clean columns -------------------
  block_out <- d$block_num
  if (renumber_blocks) {
    ## The raw platform numbers blocks continuing from any practice/
    ## instruction blocks earlier in the session (e.g. 4,5,6 instead of
    ## 1,2,3). Set renumber_blocks=TRUE to relabel them 1,2,3,... in
    ## order of appearance, matching the convention used in earlier
    ## pre-processed files in this project.
    block_out <- match(block_out, unique(block_out))
  }

  data.frame(
    PID                           = pid,
    `Task Design`                 = task_design,
    Trial                         = d$overall_idx,
    Block                         = block_out,
    `Helper Policy`               = recode_policy(d[["caring/uncaring"]]),
    `Effort Intended`             = recode_effort(d$effort_value_selected),
    `Effort Actual`               = recode_effort(d$effort_value_actual),
    `RT Effort Decision (ms)`     = d$effort_selection_rxn_time,
    `RT Satisfaction Rating (ms)` = d$rate_self_rxn_time,
    `RT Trust Rating (ms)`        = d$rate_partner_rxn_time,
    `Satisfaction Rating`         = d$rate_self_value,
    `Trust Rating`                = d$rate_partner_value,
    `Points For Helper`           = d$partner_points,
    `Points For Help-Seeker`      = d$effort_points,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


## ================================================================
##  Edit these 3 lines for each new file, then run the script
## ================================================================
if (sys.nframe() == 0) {   
  pid = 'TEST468';                      
  df <- parse_helping_task_session(
    filepath        = "helpTask2025_PARTICIPANT_SESSION_2026-07-02_16h42.46.468.csv",
    # filepath        = "helpTask2025_PARTICIPANT_SESSION_2026-06-25_14h52.52.289.csv",  # FILE NAME
    pid = pid ,
    # pid             = "LE1", # Participant ID
    task_design     = "CUC/CUU", # Task Design (Pilot was Random, Lived Exp. was CUC/CUU)
    renumber_blocks = TRUE   
  )
  write.csv(df,
            paste(pid,'.csv',sep=''),
            # "LE1.csv",
            row.names = FALSE) #Update Pts. ID
  cat("\nPreview:\n")
  print(head(df))
  cat("\nBlocks present:", paste(sort(unique(df$Block)), collapse = ", "), "\n")
}
