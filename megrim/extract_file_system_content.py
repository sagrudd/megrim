from IPython.display import Image, display, HTML, Markdown, clear_output


class extract_file_system_content:

    def __init__(self, title, object_id, path, registry=None, project_manager, aggregate=False, aggregate_key=None):
        display(Markdown("# nanopoRe: %s" % title))
        self.object_id = object_id
        self.path = path
        self.registry = registry
        self.aggregate = aggregate
        self.aggregate_key = aggregate_key
        self.project_manager = project_manager

    ##################################
    self.validate()
    self.init_components()
    self.init_interactivity()
    self.create_panel()


def validate(self):
    try:
        self.project_manager
    except NameError:
        raise RuntimeError('This workflow requires a nanopoRe project_manager object')
    else:
        if self.project_manager == None:
            raise RuntimeError('The nanopoRe project_manager has not been initialised')
    status = target_dir(self.path)


def init_components(self):
    self.sel_obj = 'select %s resource' % self.object_id
    self.web_obj = 'web accessible %s' % self.object_id
    self.obj_file = 'local %s file' % self.object_id
    self.pre_def = 'predefined %s' % self.object_id

    self.reference_id_txt = '<reference_id>'

    available_options = [self.sel_obj, self.web_obj, self.obj_file]
    if not self.registry == None:
        available_options.append(self.pre_def)

    self.reference_id = widgets.Text(value=self.reference_id_txt, disabled=False, description="name")
    self.core_object_options = widgets.Dropdown(options=available_options, description="%s" % self.object_id)
    self.object_selection_output = widgets.Output()
    self.web_obj_path = widgets.Text(value="http://", disabled=False, description="source")
    self.web_obj_button = widgets.Button(description="Upload", disabled=False, button_style="info", tooltip="upload",
                                         icon="check")
    self.file_obj_path = widgets.Text(value="/data", disabled=False, description="source")
    self.file_obj_button = widgets.Button(description="Import", disabled=False, button_style="info", tooltip="import",
                                          icon="check")
    self.existing_reference = widgets.Dropdown(options=[], description="Existing %s" % self.object_id)
    self.existing_button = widgets.Button(description="Select", disabled=False, button_style="info", tooltip="select",
                                          icon="check")
    self.button_output = widgets.Output()
    self.already_specced_output = widgets.Output()

    self.edit_setup_spec = widgets.Checkbox(value=False, description="change specifications", disabled=False)

    self.aggregate_data = widgets.Checkbox(value=False, description="aggregate data", disabled=False)


def core_object_options_change(self, change):
    if change['type'] == 'change' and change['name'] == 'value':
        with self.object_selection_output:
            clear_output()
            option = self.core_object_options.value
            # print('delta core [%s]' % option)
            if option == self.web_obj:
                if not self.registry == None:
                    display(self.reference_id)
                display(self.web_obj_path)
                display(self.web_obj_button)
                self.button_output.clear_output()
                display(self.button_output)
            elif option == self.obj_file:
                if not self.registry == None:
                    display(self.reference_id)
                display(self.file_obj_path)
                if self.aggregate:
                    display(self.aggregate_data)
                display(self.file_obj_button)
                self.button_output.clear_output()
                display(self.button_output)
            elif option == self.pre_def:
                self.existing_reference.options = self.registry.list_references()
                display(self.existing_reference)
                display(self.existing_button)
                self.button_output.clear_output()
                display(self.button_output)


def web_obj_button_clicked(self, web_obj_button):
    with self.button_output:
        option = self.web_obj_path.value
        print("web_obj_button_clicked [%s]" % option)
        dest = os.path.join(self.path, os.path.basename(option))
        if self.reference_id.value == self.reference_id_txt:
            self.reference_id.value = hashlib.md5(re.sub("^[^:]+://", "", option).encode('utf-8')).hexdigest()
        if not web_check(option):
            raise RuntimeError("specified web object [%s] cannot be found" % option)
        urllib.request.urlretrieve(option, dest, MyProgressBar())
        project_manager.set_defined_value(self.object_id, dest)
        if not self.registry == None:
            self.registry.add_item(self.reference_id.value, option, dest)


def file_obj_button_clicked(self, file_obj_button):
    with self.button_output:
        clear_output()
        option = self.file_obj_path.value
        print("file_obj_button_clicked [%s]" % option)
        if not self.aggregate_data.value:
            dest = os.path.join(self.path, os.path.basename(option))
            if self.reference_id.value == self.reference_id_txt:
                self.reference_id.value = hashlib.md5(option.encode('utf-8')).hexdigest()
            if not os.path.isfile(option):
                raise RuntimeError("specified file [%s] cannot be found" % option)
            urllib.request.urlretrieve("file://" + option, dest, MyProgressBar())
            project_manager.set_defined_value(self.object_id, dest)
            if not self.registry == None:
                self.registry.add_item(self.reference_id.value, "file://" + option, dest)
        else:
            print("aggregating data ... [key=%s]" % self.aggregate_key)
            scanner = scan_content(option, self.aggregate_key)
            fname = os.path.basename(re.sub("/*$", "", option))
            target = os.path.join(self.path, "%s.%s.gz" % (fname, self.object_id))
            scanner.aggregate(target)
            project_manager.set_defined_value(self.object_id, target)


def existing_button_clicked(self, existing_button):
    with self.button_output:
        option = self.existing_reference.value
        print("existing_button_clicked [%s]" % option)
        reference = self.registry.extract_existing(option)
        project_manager.set_defined_value(self.object_id, reference)


def edit_setup_spec_changed(self, b):
    if b['type'] == 'change':
        with self.already_specced_output:
            clear_output()
            if self.edit_setup_spec.value:
                if self.registry == None:
                    display(self.core_object_options)
                    display(self.object_selection_output)
                else:
                    self.core_object_options.value = self.pre_def
                    self.existing_reference.value = self.registry.getIdFromPath(
                        project_manager.get_defined_value(self.object_id))
                    display(self.core_object_options)
                    display(self.object_selection_output)
                    with self.object_selection_output:
                        display(self.existing_reference)
                        display(self.existing_button)
                        display(self.button_output)


def init_interactivity(self):
    self.core_object_options.observe(self.core_object_options_change)
    self.web_obj_button.on_click(self.web_obj_button_clicked)
    self.file_obj_button.on_click(self.file_obj_button_clicked)
    self.existing_button.on_click(self.existing_button_clicked)
    self.edit_setup_spec.observe(self.edit_setup_spec_changed)


def create_panel(self):
    if project_manager.has_defined_value(self.object_id):
        print("%s [%s] already selected" % (self.object_id, project_manager.get_defined_value(self.object_id)))
        display(self.edit_setup_spec)
        display(self.already_specced_output)
    else:
        print("%s needs to be created" % self.object_id)
        display(self.core_object_options)
        display(self.object_selection_output)
