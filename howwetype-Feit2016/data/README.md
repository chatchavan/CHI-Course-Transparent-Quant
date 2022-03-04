# Typing log data

The log data are based on the original log data by [Feit et al. (2016)](http://dx.doi.org/10.1145/2858036.2858233) with modifications for educational purposes by Chat Wacharamanotham.


## Changes
The following changes were made to the original raw log files:
    * Cleaning:
        * modifying characters from the `input` field that seems to be wrong. (For the full list, see `clean_map.csv`)
        * remove non-character keys from the `input` column
    * Removing redundant fields:
        * remove timestamp (kept input_time_ms)
        * remove `input` field (because it is inconsistent what this field means)
    * Clarification:
        * rename `stimulus_id` to `stimulus_index`
        * rename `current_input` to `transcribed_string`
    * Removing fields for exercise purpose:
        * remove measurement columns (so that students have to write a function to calculate WPM, etc.)
        * remove user ID and condition from the data (kept in the file name)    

## Files

**`combined_log_background.csv`** contains summarized data from the following files with the word-per-minute (WPM) calculated. If you struggle with combining raw data files together, use this file as a starting point.

**`participant_info.csv`** contains information whether the participant is a touch-typist or not

**The `log\` folder:** Each file contains typing data from one participant in one sentence condition (either "Sentence", "Mix", or "Random"). The name of the file (e.g., `User5307_T1438007839_Sentences.csv`) contains the following information:

* an anoynmized ID of study participants (e.g., `5307`)
* UNIX timestamp at the start of the condition (e.g., `1438007839`)
* and the name of the condition (e.g., `Sentences`)

Each file contains the following fields:

* `stimulus_index`: The order in which the stimulus is presented to each participant. This field is unique only within one condition of each participant. This field is *not* associated with the content of the stimulus. 
* `stimulus`: The actual sentence presented to the participants.
* `input_time_ms`: Time elapsed since the participant type the first keystroke for each sentence.
* `input_index`: The order of the keystrokes that the participants made.
* `key_symbol`: The keystrokes, including non-character keys such as Right Shift (`Shift_R`) or Backspace
* `transcribed_string`: The final text that the user produce as appear on the users' screen at the end of the trial.


## Data processing suggestions

* Use `list.files()` to list the data files for reading

* __Convert `input_time_ms` to `double` before calculation.__ Because this field is written in the log file without a decimal point, R will assume that it is an integer. When you use division in your calculation, some precision may be lost.


