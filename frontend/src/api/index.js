import axios from 'axios';

const API_BASE_URL = '/api';

// DEG Analysis API
export const degApi = {
    createJob: async (formData) => {
        const response = await axios.post(`${API_BASE_URL}/deg/jobs/`, formData, {
            headers: { 'Content-Type': 'multipart/form-data' }
        });
        return response.data;
    },

    getJob: async (jobId) => {
        const response = await axios.get(`${API_BASE_URL}/deg/jobs/${jobId}/`);
        return response.data;
    },

    getResults: async (jobId) => {
        const response = await axios.get(`${API_BASE_URL}/deg/jobs/${jobId}/results/`);
        return response.data;
    },

    listJobs: async () => {
        const response = await axios.get(`${API_BASE_URL}/deg/jobs/`);
        return response.data;
    }
};

// Annotation API
export const annotationApi = {
    mapGenes: async (genes) => {
        const response = await axios.post(`${API_BASE_URL}/annotation/gene-map/`, { genes });
        return response.data;
    },

    getKEGG: async (genes) => {
        const response = await axios.post(`${API_BASE_URL}/annotation/kegg/`, { genes });
        return response.data;
    },

    getGO: async (genes) => {
        const response = await axios.post(`${API_BASE_URL}/annotation/go/`, { genes });
        return response.data;
    },

    getDiseases: async (genes) => {
        const response = await axios.post(`${API_BASE_URL}/annotation/diseases/`, { genes });
        return response.data;
    },

    getDrugs: async (genes) => {
        const response = await axios.post(`${API_BASE_URL}/annotation/drugs/`, { genes });
        return response.data;
    }
};

// Networks API
export const networksApi = {
    getGenePathwayDisease: async (data) => {
        const response = await axios.post(`${API_BASE_URL}/networks/gene-pathway-disease/`, data);
        return response.data;
    },

    getGeneDrug: async (data) => {
        const response = await axios.post(`${API_BASE_URL}/networks/gene-drug/`, data);
        return response.data;
    },

    getTherapyNetwork: async (data) => {
        const response = await axios.post(`${API_BASE_URL}/networks/therapy/`, data);
        return response.data;
    }
};
