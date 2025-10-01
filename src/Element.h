// Element.h
#pragma once

// ******** Definition of Base Polymorphic element Class ******** //
class Element {
private:
    int id;     // Unique element identifier
public:
    // constructor & deconstructor
    Element(int elementId) : id(elementId) {}
    virtual ~Element() = default;

    // getter
    int getId() const {return id;}

    // interfaces getter via pure virtual functions
    virtual int getNnode() const = 0;
    virtual const std::vector<std::shared_ptr<Node>>& getNodes() const = 0;
    virtual double getLocalTime() const = 0;
    // interfaces setter via pure virtual functions
    virtual void setLocalTime(double t) = 0;

    // Main work: integrate / stiffness/force calculations (pure virtual)
    virtual void integrate(int kinc, int elemID, int nn, double xc[], double yc[],
                           double zc[], double disp[], double paramsR[],
                           double stf[], double qf[], double dt) = 0;
};
