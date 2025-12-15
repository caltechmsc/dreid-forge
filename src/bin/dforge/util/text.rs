pub fn wrap(text: &str, width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();

    for word in text.split_whitespace() {
        if current.is_empty() {
            current = word.to_string();
        } else if current.len() + 1 + word.len() <= width {
            current.push(' ');
            current.push_str(word);
        } else {
            lines.push(current);
            current = word.to_string();
        }
    }

    if !current.is_empty() {
        lines.push(current);
    }

    if lines.is_empty() {
        lines.push(String::new());
    }

    lines
}

pub fn truncate(s: &str, max_len: usize) -> String {
    if max_len == 0 {
        return String::new();
    }
    if max_len == 1 {
        return "…".to_string();
    }

    if s.char_indices().nth(max_len).is_none() {
        return s.to_string();
    }

    let take = max_len - 1;
    let cut = s.char_indices().nth(take).map(|(idx, _)| idx).unwrap_or(0);

    let mut out = String::with_capacity(cut + '…'.len_utf8());
    out.push_str(&s[..cut]);
    out.push('…');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_short_text() {
        let result = wrap("hello world", 20);
        assert_eq!(result, vec!["hello world"]);
    }

    #[test]
    fn wrap_long_text() {
        let result = wrap("the quick brown fox", 10);
        assert_eq!(result, vec!["the quick", "brown fox"]);
    }

    #[test]
    fn truncate_short() {
        assert_eq!(truncate("hello", 10), "hello");
    }

    #[test]
    fn truncate_exact() {
        assert_eq!(truncate("hello", 5), "hello");
    }

    #[test]
    fn truncate_long() {
        assert_eq!(truncate("hello world", 8), "hello w…");
    }

    #[test]
    fn truncate_unicode() {
        assert_eq!(truncate("日本語テスト", 4), "日本語…");
    }
}
