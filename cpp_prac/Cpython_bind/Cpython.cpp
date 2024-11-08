// g++ -o demo_module.so -shared -fPIC -I/usr/include/python3.8 Cpython.cpp

#include <Python.h>

// 简单的 C++ 函数，例如两个整数相加
static PyObject* add(PyObject* self, PyObject* args) {
    int a, b;
    // 解析 Python 传入的参数，期望两个整数
    if (!PyArg_ParseTuple(args, "ii", &a, &b)) {
        return nullptr;
    }
    int result = a + b;
    // 返回一个 Python 对象（整数）
    return PyLong_FromLong(result);
}

static PyObject* sub(PyObject* self, PyObject* args) {
    int a, b;
    // 解析 Python 传入的参数，期望两个整数
    if (!PyArg_ParseTuple(args, "ii", &a, &b)) {
        return nullptr;
    }
    int result = a - b;
    // 返回一个 Python 对象（整数）
    return PyLong_FromLong(result);
}


// 定义模块方法表
static PyMethodDef demoModuleMethods[] = {
    {"add", add, METH_VARARGS, "Add two integers"},
    {"sub", sub, METH_VARARGS, "Sub two integers"},
    {nullptr, nullptr, 0, nullptr} // 结束符
};

// 定义模块
static struct PyModuleDef demo_module = {
    PyModuleDef_HEAD_INIT,
    "demo_module",   // 模块名称
    nullptr,      // 模块文档（可选）
    -1,           // 模块状态大小（-1 表示全局模块）
    demoModuleMethods
};

// 入口函数
// 初始化函数
PyMODINIT_FUNC PyInit_demo_module(void) {
    return PyModule_Create(&demo_module);
}