"""
HTML configuration and templates for HLA-PEPCLUST reports.

This module contains HTML templates and styling for generating the results reports.
"""

# HTML template for the main page header and styling
html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HLA-PEPCLUST Results</title>
    <!-- Bootstrap CSS -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.2/dist/css/bootstrap.min.css" rel="stylesheet"
        integrity="sha384-EVSTQN3/azprG1Anm3QDgpJLIm9Nao0Yz1ztcQTwFspd3yD65VohhpuuCOmLASjC" crossorigin="anonymous">
    <!-- Bootstrap Icons -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.1/font/bootstrap-icons.css">
    <!-- DataTables CSS -->
    <link href="https://cdn.datatables.net/2.2.2/css/dataTables.bootstrap5.css" rel="stylesheet"/>
    <!-- Vega & Vega-Lite -->
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    
    <style>
        /* Custom styles */
        .no-img {
            display: block;
            width: 100%;
            height: 260px;
            background-color: #f8f9fa;
            border: 1px dashed #dee2e6;
            text-align: center;
            line-height: 260px;
            color: #6c757d;
        }
        .card {
            margin-bottom: 20px;
            box-shadow: 0 0.125rem 0.25rem rgba(0,0,0,0.075);
        }
        .card-img {
            height: 260px;
            object-fit: contain;
            padding: 10px;
        }
        h5 {
            padding: 10px;
            margin: 0;
            background-color: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
        }
        .row {
            margin-bottom: 20px;
        }
    </style>
</head>
<body>
"""

# HTML template for the body content start
body_start = """
<header>
    <div class="navbar navbar-dark bg-dark shadow-sm">
        <div class="container">
            <a href="#" class="navbar-brand d-flex align-items-center">
                <i class="bi bi-diagram-3 me-2"></i>
                <strong>HLA-PEPCLUST Results</strong>
            </a>
        </div>
    </div>
</header>

<main>
    <section class="py-5 container">
        <div class="row">
            <div class="col-lg-12 mx-auto">
                <h1 class="fw-light">Cluster Search Results</h1>
                <p class="lead text-muted">
                    This report shows the correlation between input clusters and reference HLA/MHC motifs.
                </p>
            </div>
        </div>
        
        <div class="container">
            <div class="row">
"""

# Functions to generate HTML components will be imported from cluster_search.py