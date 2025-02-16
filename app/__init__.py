# app/__init__.py
from flask import Flask
from .celery_tasks import make_celery


def create_app():
    app = Flask(__name__)

    # Load configuration from the Config class or directly
    app.config.from_object("config.Config")

    # Initialize Celery
    celery = make_celery(app)

    # Register Blueprints and other routes
    from .routes import bp

    app.register_blueprint(bp)

    return app
