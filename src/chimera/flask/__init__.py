from flask import Flask

from chimera import __version__, config


def create_app(debug=False):

    app = Flask('chimera', static_folder='flask/static', template_folder='flask/templates')
    app.config.from_object(config.flask)

    if debug:
        from werkzeug.debug import DebuggedApplication
        app.wsgi_app = DebuggedApplication(app.wsgi_app, True)

    import chimera.flask.blueprints.web

    app.register_blueprint(chimera.flask.blueprints.web.bp, url_prefix='/')

    app.add_template_global(lambda: __version__, name='app_version')
    return app
