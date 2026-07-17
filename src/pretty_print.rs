// stdout styles
const BOLD: &str = "\x1b[0;1m";
const ITALIC: &str = "\x1b[0;3m";
const UNDERL: &str = "\x1b[0;4m";
const RED: &str = "\x1b[91m";
const GREEN: &str = "\x1b[92m";
const BLUE: &str = "\x1b[94m";
const RED_BG: &str = "\x1b[41m";
const GREEN_BG: &str = "\x1b[42m";
const BLUE_BG: &str = "\x1b[44m";
const NORMAL: &str = "\x1b[0m";

pub fn bold(s: &str) -> String {
    format!("{BOLD}{s}{NORMAL}")
}

pub fn italic(s: &str) -> String {
    format!("{ITALIC}{s}{NORMAL}")
}

pub fn underl(s: &str) -> String {
    format!("{UNDERL}{s}{NORMAL}")
}

pub fn red(s: &str) -> String {
    format!("{RED}{s}{NORMAL}")
}

pub fn blue(s: &str) -> String {
    format!("{BLUE}{s}{NORMAL}")
}

pub fn green(s: &str) -> String {
    format!("{GREEN}{s}{NORMAL}")
}

pub fn red_bg(s: &str) -> String {
    format!("{RED_BG}{s}{NORMAL}")
}

pub fn blue_bg(s: &str) -> String {
    format!("{BLUE_BG}{s}{NORMAL}")
}

pub fn green_bg(s: &str) -> String {
    format!("{GREEN_BG}{s}{NORMAL}")
}

