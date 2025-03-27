# Mkdocs documentation

This file describes the Mkdocs package we use inside Gyselalib to document our code.

## Packages for installation

The packages required to make the docs are listed in the [requirements.txt](./requirements.txt) file. Simply run:

```bash
pip install -r requirements.txt
```

Adjustments to MkDocs can be done in the `mkdocs.yml`.

To run the documentation locally follow the instructions in [Add Documentation](./development/Adding_docs.md).

## MkDocs Navigation and Configuration

The configuration of MkDocs, including various options and additional packages, is managed through the `mkdocs.yml` file. This file allows you to customize the structure and appearance of your documentation.

### Modifying the Navigation Panel

To add a new section to the navigation panel of your documentation site, update the `nav:` section in the `mkdocs.yml` file. This section defines the structure of your documentation, specifying the order and hierarchy of pages displayed in the sidebar.

### Configuring the Doxygen API

The documentation includes an API reference generated with Doxygen, you can configure it using [`mkdoxy`](https://mkdoxy.kubaandrysek.cz/), an additional package designed to integrate Doxygen into MkDocs.

## Useful pages used for setting up MkDocs

- find a nice example project here: [MkDocs Project](https://example-mkdocs-basic.readthedocs.io/en/latest/#example-project-usage)
- mkdocs-include-dir-to-nav ... Package to do include subdirectories [Doc](https://github.com/mysiki/mkdocs_include_dir_to_nav)
