const express = require('express');
const cors = require('cors');
const bodyParser = require('body-parser');
const fs = require('fs');
const path = require('path');

const app = express();
const PORT = 3000;

// Middleware
app.use(cors());
app.use(bodyParser.json({ limit: '50mb' }));
app.use(bodyParser.urlencoded({ extended: true, limit: '50mb' }));

// Хранилище данных (в реальном приложении используйте БД)
let analysisResults = [];
let nextId = 1;

// Маршруты API

// Получить все результаты анализа
app.get('/api/analysis', (req, res) => {
    res.json({
        success: true,
        data: analysisResults,
        count: analysisResults.length
    });
});

// Получить конкретный результат по ID
app.get('/api/analysis/:id', (req, res) => {
    const id = parseInt(req.params.id);
    const result = analysisResults.find(item => item.id === id);
    
    if (!result) {
        return res.status(404).json({
            success: false,
            message: 'Analysis result not found'
        });
    }
    
    res.json({
        success: true,
        data: result
    });
});

// Отправить новые данные анализа от сайта
app.post('/api/analysis', (req, res) => {
    try {
        const analysisData = req.body;
        
        // Валидация данных
        if (!analysisData.mesh || !analysisData.stats) {
            return res.status(400).json({
                success: false,
                message: 'Invalid data format: mesh and stats are required'
            });
        }
        
        const newResult = {
            id: nextId++,
            timestamp: new Date().toISOString(),
            ...analysisData
        };
        
        analysisResults.push(newResult);
        
        // Сохраняем в файл для C++ приложения
        saveToFile(newResult);
        
        res.json({
            success: true,
            message: 'Analysis data received successfully',
            id: newResult.id
        });
        
    } catch (error) {
        console.error('Error processing analysis data:', error);
        res.status(500).json({
            success: false,
            message: 'Internal server error'
        });
    }
});

// Экспорт данных в формате для C++
app.get('/api/export/cpp/:id', (req, res) => {
    const id = parseInt(req.params.id);
    const result = analysisResults.find(item => item.id === id);
    
    if (!result) {
        return res.status(404).json({
            success: false,
            message: 'Analysis result not found'
        });
    }
    
    // Форматируем данные для C++ приложения
    const cppData = formatForCpp(result);
    
    res.json({
        success: true,
        data: cppData
    });
});

// Удалить результат анализа
app.delete('/api/analysis/:id', (req, res) => {
    const id = parseInt(req.params.id);
    const index = analysisResults.findIndex(item => item.id === id);
    
    if (index === -1) {
        return res.status(404).json({
            success: false,
            message: 'Analysis result not found'
        });
    }
    
    analysisResults.splice(index, 1);
    
    res.json({
        success: true,
        message: 'Analysis result deleted successfully'
    });
});

// Вспомогательные функции
function saveToFile(data) {
    const filename = `analysis_${data.id}_${Date.now()}.json`;
    const filepath = path.join(__dirname, 'exports', filename);
    
    // Создаем папку exports если её нет
    const exportsDir = path.join(__dirname, 'exports');
    if (!fs.existsSync(exportsDir)) {
        fs.mkdirSync(exportsDir);
    }
    
    fs.writeFileSync(filepath, JSON.stringify(data, null, 2));
    console.log(`Data saved to ${filename}`);
}

function formatForCpp(data) {
    return {
        analysis_id: data.id,
        timestamp: data.timestamp,
        mesh_info: {
            node_count: data.stats.node_count,
            element_count: data.stats.element_count,
            element_type: data.stats.element_type
        },
        statistics: data.stats,
        bounds: data.bounds,
        metadata: data.metadata || {}
    };
}

// Запуск сервера
app.listen(PORT, () => {
    console.log(`REST API server running on http://localhost:${PORT}`);
    console.log(`C++ exports directory: ${path.join(__dirname, 'exports')}`);
});