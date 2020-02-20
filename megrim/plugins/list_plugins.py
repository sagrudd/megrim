from megrim.environment import MegrimPlugin


class ListPlugins(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "list"

    def execute(self, args):
        print(21)

    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="list stuff help", parents=[parent_parser])

