import threading
import unittest

from inchworm_mvp import expand_forward


class ExpandForwardTests(unittest.TestCase):
    def test_does_not_loop_forever_on_repeated_kmers(self):
        kmers_count = {"AAA": 3}
        result_container: dict[str, str] = {}

        def run_function() -> None:
            result_container["value"] = expand_forward(kmers_count, "AAA", 3)

        thread = threading.Thread(target=run_function, daemon=True)
        thread.start()
        thread.join(0.5)

        self.assertFalse(thread.is_alive(), "expand_forward should stop instead of looping forever")
        self.assertEqual(result_container.get("value", ""), "AAAA")


if __name__ == "__main__":
    unittest.main()
