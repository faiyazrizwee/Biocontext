from django.urls import path
from .views import (
    GenePathwayDiseaseNetworkView,
    GeneDrugNetworkView,
    FullTherapyNetworkView
)

urlpatterns = [
    path('gene-pathway-disease/', GenePathwayDiseaseNetworkView.as_view(), name='gene-pathway-disease-network'),
    path('gene-drug/', GeneDrugNetworkView.as_view(), name='gene-drug-network'),
    path('therapy/', FullTherapyNetworkView.as_view(), name='full-therapy-network'),
]
