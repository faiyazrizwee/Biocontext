from django.urls import path
from .views import (
    GeneMapView,
    KEGGEnrichmentView,
    GOEnrichmentView,
    DiseaseAssociationView,
    DrugAssociationView
)

urlpatterns = [
    path('gene-map/', GeneMapView.as_view(), name='gene-map'),
    path('kegg/', KEGGEnrichmentView.as_view(), name='kegg-enrichment'),
    path('go/', GOEnrichmentView.as_view(), name='go-enrichment'),
    path('diseases/', DiseaseAssociationView.as_view(), name='disease-association'),
    path('drugs/', DrugAssociationView.as_view(), name='drug-association'),
]
