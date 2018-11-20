HEADERS += \
    test_function.h \
    hybrid_method.h \
    types.h \
    test_function_collection.h \
    tf_registry.h

SOURCES += \
    test_function.cpp \
    main.cpp \
    hybrid_method.cpp \
    test_functions/tf_diagonal.cpp \
    test_functions/tf_extended.cpp \
    test_functions/tf_full_hessian.cpp \
    test_functions/tf_quadratic.cpp \
    test_functions/tf_cute.cpp \
    test_functions/tf_generalized_.cpp \
    tf_registry.cpp

QT -= gui core
QMAKE_CXX += -std=c++14 -O3
QMAKE_LFLAGS += -pipe -Wall -fPIC -v
LIBS += -larmadillo   -L/usr/local/lib/libopenblas.a -llapack -fopenmp
