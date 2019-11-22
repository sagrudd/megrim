import ipywidgets as widgets
from IPython.display import display, clear_output

from .project_manager import project_manager


class project_selector:

    def __init__(self, epi2me_home, projects, project_type, project_registry):

        self.epi2me_home = epi2me_home
        self.projects = projects
        self.project_type = project_type
        self.pr = project_registry
        self.project_manager = None

        self.project_panel = widgets.Output()
        self.project_selection_panel = widgets.Output()

        self.project_item_select = "Project Select"
        self.project_item_new = "New Project"
        self.project_item_existing = "Existing Project"

        self.project_definition = widgets.Dropdown(
            options=[self.project_item_select, self.project_item_new, self.project_item_existing],
            description="choose project")
        self.existing_projects = widgets.Dropdown(options=[], description="choose project")

        self.project_name_txt = "<specify project name>"
        self.project_name = widgets.Text(value=self.project_name_txt, disabled=False)
        self.add_project_button = widgets.Button(description="Add project", disabled=False, button_style="info",
                                                 tooltip="Click to add", icon="check")
        self.existing_project_button = widgets.Button(description="Use project", disabled=False, button_style="info",
                                                      tooltip="Click to use", icon="check")

        self.project_definition.observe(self.project_definition_change)
        self.add_project_button.on_click(self.add_project_button_clicked)
        self.existing_project_button.on_click(self.existing_project_button_clicked)

        with self.project_panel:
            display(self.project_definition)
            display(self.project_selection_panel)

    def add_project_button_clicked(self, button):
        with self.project_selection_panel:
            option = removeDisallowedFilenameChars(self.project_name.value)
            print("add button clicked [%s] ..." % option)
            if option == self.project_name_txt:
                raise RuntimeError('please specify a project name for this [%s] study' % self.project_type)
            if not target_dir(os.path.join(self.epi2me_home, self.projects, option)):
                raise RuntimeError(
                    'project [%s] already exists ... aborting ' % os.path.join(epi2me_home, projects, option))
            self.pr.add_project(option, os.path.join(self.epi2me_home, self.projects, option), self.project_type)
            project_path = self.pr.get_project_path(option)
            self.project_manager = project_manager(project_path, self.project_type, option)
            print(project_path)

    def existing_project_button_clicked(self, button):
        with self.project_selection_panel:
            option = self.existing_projects.value
            print("use button clicked [%s] ..." % option)
            project_path = self.pr.get_project_path(option)
            self.project_manager = project_manager(project_path, self.project_type, option)
            print(project_path)

    def project_definition_change(self, change):
        if change['type'] == 'change' and change['name'] == 'value':
            with self.project_selection_panel:
                clear_output()
                option = self.project_definition.value

                if option == self.project_item_new:
                    print("option == [%s]" % option)
                    display(self.project_name)
                    display(self.add_project_button)
                elif option == self.project_item_existing:
                    print("option == [%s]" % option)
                    self.existing_projects.options = self.pr.list_projects(self.project_type)
                    display(self.existing_projects)
                    display(self.existing_project_button)

    def get_project_manager(self):
        return self.project_manager

    def display(self):
        display(self.project_panel)
