import unittest

from pathway_mvp import parse_kegg_tsv


class ParseKeggTsvTests(unittest.TestCase):
    def test_parses_one_pair_per_line(self):
        text = "hsa:6772\tpath:hsa04630\nhsa:6772\tpath:hsa04658\n"
        self.assertEqual(
            parse_kegg_tsv(text),
            [("hsa:6772", "path:hsa04630"), ("hsa:6772", "path:hsa04658")],
        )

    def test_ignores_blank_lines(self):
        text = "hsa:6772\tpath:hsa04630\n\nhsa:6772\tpath:hsa04658\n"
        self.assertEqual(len(parse_kegg_tsv(text)), 2)

    def test_returns_empty_list_for_empty_response(self):
        self.assertEqual(parse_kegg_tsv(""), [])

    def test_value_containing_a_tab_is_kept_whole(self):
        # KEGG's `list` response values (pathway names) can contain further
        # tabs/formatting - only the first tab separates key from value.
        text = "path:hsa04630\tJAK-STAT signaling pathway\t- Homo sapiens (human)\n"
        self.assertEqual(
            parse_kegg_tsv(text),
            [("path:hsa04630", "JAK-STAT signaling pathway\t- Homo sapiens (human)")],
        )


if __name__ == "__main__":
    unittest.main()
