/* lsqcpp.h
 *
 * Author: Fabian Meyer
 * Created On: 22 Jul 2019
 * License: MIT
 */

#ifndef LSQCPP_LSQCPP_H_
#define LSQCPP_LSQCPP_H_

#include <gsEigen/Geometry>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <functional>

namespace lsq
{
    typedef gsEigen::MatrixXd::Index Index;

    /** Functor to compute forward differences.
      * Computes the gradient of the objective f(x) as follows:
      *
      * grad(x) = (f(x + eps) - f(x)) / eps
      *
      * The computation requires len(x) evaluations of the objective.
      */
    template<typename Scalar>
    class ForwardDifferences
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &)> ErrorFunction;
    private:
        Scalar eps_;
        int threads_;
        ErrorFunction objective_;
    public:
        ForwardDifferences()
            : ForwardDifferences(
                std::sqrt(std::numeric_limits<Scalar>::epsilon()))
        { }

        ForwardDifferences(const Scalar eps)
            : eps_(eps), threads_(1), objective_()
        { }

        void setNumericalEpsilon(const Scalar eps)
        {
            eps_ = eps;
        }

        void setThreads(const int threads)
        {
            threads_ = threads;
        }

        void setErrorFunction(const ErrorFunction &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Vector &fval,
            Matrix &jacobian)
        {
            assert(objective_);

            jacobian.resize(fval.size(), xval.size());
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < xval.size(); ++i)
            {
                Vector fvalN;
                Vector xvalN = xval;
                xvalN(i) += eps_;
                objective_(xvalN, fvalN);

                jacobian.col(i) = (fvalN - fval) / eps_;
            }
        }
    };

    /** Functor to compute backward differences.
      * Computes the gradient of the objective f(x) as follows:
      *
      * grad(x) = (f(x) - f(x - eps)) / eps
      *
      * The computation requires len(x) evaluations of the objective.
      */
    template<typename Scalar>
    class BackwardDifferences
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &)> ErrorFunction;
    private:
        Scalar eps_;
        int threads_;
        ErrorFunction objective_;
    public:
        BackwardDifferences()
            : BackwardDifferences(
                std::sqrt(std::numeric_limits<Scalar>::epsilon()))
        { }

        BackwardDifferences(const Scalar eps)
            : eps_(eps), threads_(1), objective_()
        { }

        void setNumericalEpsilon(const Scalar eps)
        {
            eps_ = eps;
        }

        void setThreads(const int threads)
        {
            threads_ = threads;
        }

        void setErrorFunction(const ErrorFunction &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Vector &fval,
            Matrix &jacobian)
        {
            assert(objective_);

            jacobian.resize(fval.size(), xval.size());
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < xval.size(); ++i)
            {
                Vector fvalN;
                Vector xvalN = xval;
                xvalN(i) -= eps_;
                objective_(xvalN, fvalN);
                jacobian.col(i) = (fval - fvalN) / eps_;
            }
        }
    };

    /** Functor to compute central differences.
      * Computes the gradient of the objective f(x) as follows:
      *
      * grad(x) = (f(x + 0.5 eps) - f(x - 0.5 eps)) / eps
      *
      * The computation requires 2 * len(x) evaluations of the objective.
      */
    template<typename Scalar>
    struct CentralDifferences
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &)> ErrorFunction;
    private:
        Scalar eps_;
        int threads_;
        ErrorFunction objective_;
    public:
        CentralDifferences()
            : CentralDifferences(
                std::sqrt(std::numeric_limits<Scalar>::epsilon()))
        { }

        CentralDifferences(const Scalar eps)
            : eps_(eps), threads_(1), objective_()
        { }

        void setNumericalEpsilon(const Scalar eps)
        {
            eps_ = eps;
        }

        void setThreads(const int threads)
        {
            threads_ = threads;
        }

        void setErrorFunction(const ErrorFunction &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Vector &fval,
            Matrix &jacobian)
        {
            assert(objective_);

            std::vector<Vector> fvalN(xval.size() * 2);
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < static_cast<Index>(fvalN.size()); ++i)
            {
                Index idx = i / 2;
                Vector xvalN = xval;
                if(i % 2 == 0)
                    xvalN(idx) += eps_ / 2;
                else
                    xvalN(idx) -= eps_ / 2;

                objective_(xvalN, fvalN[i]);
            }

            jacobian.resize(fval.size(), xval.size());
            for(Index i = 0; i < xval.size(); ++i)
                jacobian.col(i) = (fvalN[i * 2] - fvalN[i * 2 + 1]) / eps_;
        }
    };

    /** Dummy callback functor, which does nothing. */
    template<typename Scalar>
    struct NoCallback
    {
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;

        bool operator()(const Index,
            const Vector &,
            const Vector &,
            const Matrix &,
            const Vector &,
            const Vector &) const
        {
            return true;
        }
    };

    /** Step size functor, which returns a constant step size. */
    template<typename Scalar>
    class ConstantStepSize
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &, Matrix &)> ErrorFunction;
    private:
        Scalar stepSize_;
    public:

        ConstantStepSize()
            : ConstantStepSize(static_cast<Scalar>(1.0))
        {

        }

        ConstantStepSize(const Scalar stepSize)
            : stepSize_(stepSize)
        {

        }

        /** Set the step size returned by this functor.
          * @param stepSize step size returned by functor */
        void setStepSize(const Scalar stepSize)
        {
            stepSize_ = stepSize;
        }

        void setErrorFunction(const ErrorFunction &)
        { }

        Scalar operator()(const Vector &,
            const Vector &,
            const Matrix &,
            const Vector &,
            const Vector &)
        {
            return stepSize_;
        }
    };

    /** Step size functor to compute Barzilai-Borwein (BB) steps.
      * The functor can either compute the direct or inverse BB step.
      * The steps are computed as follows:
      *
      * s_k = x_k - x_k-1         k >= 1
      * y_k = step_k - step_k-1   k >= 1
      * Direct:  stepSize = (s_k^T * s_k) / (y_k^T * s_k)
      * Inverse: stepSize = (y_k^T * s_k) / (y_k^T * y_k)
      *
      * The very first step is computed as a constant. */
    template<typename Scalar>
    class BarzilaiBorwein
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &, Matrix &)> ErrorFunction;

        enum class Method
        {
            Direct,
            Inverse
        };
    private:
        Vector lastXval_;
        Vector lastStep_;
        Method method_;
        Scalar constStep_;

        Scalar constantStep() const
        {
            return constStep_;
        }

        Scalar directStep(const Vector &xval,
            const Vector &step)
        {
            auto sk = xval - lastXval_;
            auto yk = step - lastStep_;
            Scalar num = sk.dot(sk);
            Scalar denom = sk.dot(yk);

            if(denom == 0)
                return 1;
            else
                return num / denom;
        }

        Scalar inverseStep(const Vector &xval,
            const Vector &step)
        {
            auto sk = xval - lastXval_;
            auto yk = step - lastStep_;
            Scalar num = sk.dot(yk);
            Scalar denom = yk.dot(yk);

            if(denom == 0)
                return 1;
            else
                return num / denom;
        }
    public:
        BarzilaiBorwein()
            : BarzilaiBorwein(Method::Direct, static_cast<Scalar>(1e-2))
        { }

        BarzilaiBorwein(const Method method, const Scalar constStep)
            : lastXval_(), lastStep_(), method_(method),
            constStep_(constStep)
        { }

        void setErrorFunction(const ErrorFunction &)
        { }

        void setMethod(const Method method)
        {
            method_ = method;
        }

        void setConstStepSize(const Scalar stepSize)
        {
            constStep_ = stepSize;
        }

        Scalar operator()(const Vector &xval,
            const Vector &,
            const Matrix &,
            const Vector &,
            const Vector &step)
        {
            Scalar stepSize = 0;
            if(lastXval_.size() == 0)
            {
                stepSize = (1 / step.norm()) * constStep_;
            }
            else
            {
                switch(method_)
                {
                case Method::Direct:
                    stepSize = directStep(xval, step);
                    break;
                case Method::Inverse:
                    stepSize = inverseStep(xval, step);
                    break;
                default:
                    assert(false);
                    break;
                }
            }

            lastStep_ = step;
            lastXval_ = xval;

            return stepSize;
        }
    };

    /** Step size functor to perform Armijo Linesearch with backtracking.
      * The functor iteratively decreases the step size until the following
      * conditions are met:
      *
      * Armijo: f(x - stepSize * grad(x)) <= f(x) - c1 * stepSize * grad(x)^T * grad(x)
      *
      * If the condition does not hold the step size is decreased:
      *
      * stepSize = decrease * stepSize
      */
    template<typename Scalar>
    class ArmijoBacktracking
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &, Matrix &)> ErrorFunction;
    private:
        Scalar decrease_;
        Scalar c1_;
        Scalar minStep_;
        Scalar maxStep_;
        Index maxIt_;
        ErrorFunction objective_;

    public:
        ArmijoBacktracking()
            : ArmijoBacktracking(
                static_cast<Scalar>(0.8),
                static_cast<Scalar>(1e-4),
                static_cast<Scalar>(1e-10),
                static_cast<Scalar>(1.0),
                0)
        { }

        ArmijoBacktracking(const Scalar decrease,
            const Scalar c1,
            const Scalar minStep,
            const Scalar maxStep,
            const Index iterations)
            : decrease_(decrease), c1_(c1), minStep_(minStep),
            maxStep_(maxStep), maxIt_(iterations), objective_()
        { }

        /** Set the decreasing factor for backtracking.
          * Assure that decrease in (0, 1).
          * @param decrease decreasing factor */
        void setBacktrackingDecrease(const Scalar decrease)
        {
            assert(decrease > static_cast<Scalar>(0));
            assert(decrease < static_cast<Scalar>(1));
            decrease_ = decrease;
        }

        /** Set the relaxation constant for the Armijo condition (see class description).
          * Typically c1 is chosen to be quite small, e.g. 1e-4.
          * Assure that c1 in (0, 0.5).
          * @param c1 armijo constant */
        void setArmijoConstant(const Scalar c1)
        {
            assert(c1 > static_cast<Scalar>(0));
            assert(c1 < static_cast<Scalar>(0.5));
            c1_ = c1;
        }

        /** Set the bounds for the step size during linesearch.
          * The final step size is guaranteed to be in [minStep, maxStep].
          * The
          * @param minStep minimum step size
          * @param maxStep maximum step size */
        void setStepBounds(const Scalar minStep, const Scalar maxStep)
        {
            assert(minStep < maxStep);
            minStep_ = minStep;
            maxStep_ = maxStep;
        }

        /** Set the maximum number of iterations.
          * Set to 0 or negative for infinite iterations.
          * @param iterations maximum number of iterations */
        void setMaxIterations(const Index iterations)
        {
            maxIt_ = iterations;
        }

        void setErrorFunction(const ErrorFunction &objective)
        {
            objective_ = objective;
        }

        Scalar operator()(const Vector &xval,
            const Vector &fval,
            const Matrix &,
            const Vector &gradient,
            const Vector &step)
        {
            assert(objective_);

            Scalar stepSize = maxStep_ / decrease_;
            Matrix jacobianN;
            Vector gradientN;
            Vector xvalN;
            Vector fvalN;

            Scalar error = static_cast<Scalar>(0.5) * fval.squaredNorm();
            Scalar stepGrad = gradient.dot(step);
            bool armijoCondition = false;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                stepSize * decrease_ >= minStep_ &&
                !armijoCondition)
            {
                stepSize = decrease_ * stepSize;
                xvalN = xval - stepSize * step;
                objective_(xvalN, fvalN, jacobianN);
                Scalar errorN = static_cast<Scalar>(0.5) * fvalN.squaredNorm();
                gradientN = jacobianN.transpose() * fvalN;

                armijoCondition = errorN <= error + c1_ * stepSize * stepGrad;

                ++iterations;
            }

            return stepSize;
        }
    };

    /** Step size functor to perform Wolfe Linesearch with backtracking.
      * The functor iteratively decreases the step size until the following
      * conditions are met:
      *
      * Armijo: f(x - stepSize * grad(x)) <= f(x) - c1 * stepSize * grad(x)^T * grad(x)
      * Wolfe: grad(x)^T grad(x - stepSize * grad(x)) <= c2 * grad(x)^T * grad(x)
      *
      * If either condition does not hold the step size is decreased:
      *
      * stepSize = decrease * stepSize
      */
    template<typename Scalar>
    class WolfeBacktracking
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef std::function<void(const Vector &, Vector &, Matrix &)> ErrorFunction;
    private:
        Scalar decrease_;
        Scalar c1_;
        Scalar c2_;
        Scalar minStep_;
        Scalar maxStep_;
        Index maxIt_;
        ErrorFunction objective_;

    public:
        WolfeBacktracking()
            : WolfeBacktracking(
                static_cast<Scalar>(0.8),
                static_cast<Scalar>(1e-4),
                static_cast<Scalar>(0.9),
                static_cast<Scalar>(1e-10),
                static_cast<Scalar>(1.0),
                0)
        { }

        WolfeBacktracking(const Scalar decrease,
            const Scalar c1,
            const Scalar c2,
            const Scalar minStep,
            const Scalar maxStep,
            const Index iterations)
            : decrease_(decrease), c1_(c1), c2_(c2), minStep_(minStep),
            maxStep_(maxStep), maxIt_(iterations), objective_()
        { }

        /** Set the decreasing factor for backtracking.
          * Assure that decrease in (0, 1).
          * @param decrease decreasing factor */
        void setBacktrackingDecrease(const Scalar decrease)
        {
            decrease_ = decrease;
        }

        /** Set the wolfe constants for Armijo and Wolfe condition (see class
          * description).
          * Assure that c1 < c2 < 1 and c1 in (0, 0.5).
          * Typically c1 is chosen to be quite small, e.g. 1e-4.
          * @param c1 armijo constant
          * @param c2 wolfe constant */
        void setWolfeConstants(const Scalar c1, const Scalar c2)
        {
            assert(c1 > static_cast<Scalar>(0));
            assert(c1 < static_cast<Scalar>(0.5));
            assert(c1 < c2);
            assert(c2 < static_cast<Scalar>(1));
            c1_ = c1;
            c2_ = c2;
        }

        /** Set the bounds for the step size during linesearch.
          * The final step size is guaranteed to be in [minStep, maxStep].
          * @param minStep minimum step size
          * @param maxStep maximum step size */
        void setStepBounds(const Scalar minStep, const Scalar maxStep)
        {
            assert(minStep < maxStep);
            minStep_ = minStep;
            maxStep_ = maxStep;
        }

        /** Set the maximum number of iterations.
          * Set to 0 or negative for infinite iterations.
          * @param iterations maximum number of iterations */
        void setMaxIterations(const Index iterations)
        {
            maxIt_ = iterations;
        }

        void setErrorFunction(const ErrorFunction &objective)
        {
            objective_ = objective;
        }

        Scalar operator()(const Vector &xval,
            const Vector &fval,
            const Matrix &,
            const Vector &gradient,
            const Vector &step)
        {
            assert(objective_);

            Scalar stepSize = maxStep_ / decrease_;
            Matrix jacobianN;
            Vector gradientN;
            Vector xvalN;
            Vector fvalN;

            Scalar error = fval.squaredNorm() / 2;
            Scalar stepGrad = gradient.dot(step);
            bool armijoCondition = false;
            bool wolfeCondition = false;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                stepSize * decrease_ >= minStep_ &&
                !(armijoCondition && wolfeCondition))
            {
                stepSize = decrease_ * stepSize;
                xvalN = xval - stepSize * step;
                objective_(xvalN, fvalN, jacobianN);
                Scalar errorN = fvalN.squaredNorm() / 2;
                gradientN = jacobianN.transpose() * fvalN;

                armijoCondition = errorN <= error + c1_ * stepSize * stepGrad;
                wolfeCondition = gradientN.dot(step) >= c2_ * stepGrad;

                ++iterations;
            }

            return stepSize;
        }
    };

    template<typename Scalar>
    struct DenseSVDSolver
    {
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef gsEigen::JacobiSVD<Matrix, gsEigen::FullPivHouseholderQRPreconditioner>
            Solver;

        void operator()(const Matrix &A, const Vector &b, Vector &result) const
        {
            Solver solver(A, gsEigen::ComputeFullU | gsEigen::ComputeFullV);
            result = solver.solve(b);
        }
    };

    template<typename Scalar>
    struct DenseCholeskySolver
    {
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
        typedef gsEigen::LDLT<Matrix> Solver;

        void operator()(const Matrix &A, const Vector &b, Vector &result) const
        {
            Solver decomp;
            decomp.compute(A);

            if(!decomp.isPositive())
                throw std::runtime_error("DenseCholeskySolver: matrix is not positive semi-definite");

            result = decomp.solve(b);
        }
    };

    /** Base class for least squares algorithms.
      * It implements the whole optimization strategy except the step
      * calculation. Cannot be instantiated. */
    template<typename Scalar,
        typename ErrorFunction,
        typename StepSize,
        typename Callback,
        typename FiniteDifferences>
    class LeastSquaresAlgorithm
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
    protected:
        ErrorFunction errorFunction_;
        StepSize stepSize_;
        Callback callback_;
        FiniteDifferences finiteDifferences_;

        Index maxIt_;
        Scalar minStepLen_;
        Scalar minGradLen_;
        Scalar minError_;
        Index verbosity_;

        void evaluateErrorFunction(const Vector &xval, Vector &fval, Matrix &jacobian)
        {
            jacobian.resize(0, 0);
            errorFunction_(xval, fval, jacobian);
            if(jacobian.size() == 0)
                finiteDifferences_(xval, fval, jacobian);
        }

        std::string vector2str(const Vector &vec) const
        {
            std::stringstream ss1;
            ss1 << std::fixed << std::showpoint << std::setprecision(6);
            std::stringstream ss2;
            ss2 << '[';
            for(Index i = 0; i < vec.size(); ++i)
            {
                ss1 << vec(i);
                ss2 << std::setfill(' ') << std::setw(10) << ss1.str();
                if(i != vec.size() - 1)
                    ss2 << ' ';
                ss1.str("");
            }
            ss2 << ']';

            return ss2.str();
        }

    public:
        struct Result
        {
            Vector xval;
            Vector fval;
            Scalar error;
            Index iterations;
            bool converged;
        };

        LeastSquaresAlgorithm()
            : errorFunction_(),
            stepSize_(),
            callback_(),
            finiteDifferences_(),
            maxIt_(0),
            minStepLen_(static_cast<Scalar>(1e-9)),
            minGradLen_(static_cast<Scalar>(1e-9)),
            minError_(static_cast<Scalar>(0)),
            verbosity_(0)
        { }

        virtual ~LeastSquaresAlgorithm()
        { }

        /** Set the number of threads used to compute gradients.
          * This only works if OpenMP is enabled.
          * Set to 0 to allow automatic detection of thread number.
          * @param threads number of threads to be used */
        void setThreads(const Index threads)
        {
            finiteDifferences_.setThreads(threads);
        }

        /** Set the difference for gradient estimation with finite differences.
          * @param eps numerical epsilon */
        void setNumericalEpsilon(const Scalar eps)
        {
            finiteDifferences_.setNumericalEpsilon(eps);
        }

        /** Sets the instance values of the custom error function.
          * Should be used if the error function requires custom data parameters.
          * @param errorFunction instance that should be copied */
        void setErrorFunction(const ErrorFunction &errorFunction)
        {
            errorFunction_ = errorFunction;
        }

        void setCallback(const Callback &callback)
        {
            callback_ = callback;
        }

        /** Sets the instance values of the step size functor.
          * @param stepSize instance that should be copied */
        void setStepSize(const StepSize &stepSize)
        {
            stepSize_ = stepSize;
        }

        /** Set the maximum number of iterations.
          * Set to 0 or negative for infinite iterations.
          * @param iterations maximum number of iterations */
        void setMaxIterations(const Index iterations)
        {
            maxIt_ = iterations;
        }

        /** Set the minimum step length between two iterations.
          * If the step length falls below this value, the optimizer stops.
          * @param steplen minimum step length */
        void setMinStepLength(const Scalar steplen)
        {
            minStepLen_ = steplen;
        }

        /** Set the minimum gradient length.
          * If the gradient length falls below this value, the optimizer stops.
          * @param gradlen minimum gradient length */
        void setMinGradientLength(const Scalar gradlen)
        {
            minGradLen_ = gradlen;
        }

        /** Set the minimum squared error.
          * If the error falls below this value, the optimizer stops.
          * @param error minimum error */
        void setMinError(const Scalar error)
        {
            minError_ = error;
        }

        /** Set the level of verbosity to print status information after each
          * iteration.
          * Set to 0 to deacticate any output.
          * @param verbosity level of verbosity */
        void setVerbosity(const Index verbosity)
        {
            verbosity_ = verbosity;
        }

        Result minimize(const Vector &initialGuess)
        {
            finiteDifferences_.setErrorFunction(
                [this](const Vector &xval, Vector &fval)
                { Matrix tmp; this->errorFunction_(xval, fval, tmp); });
            stepSize_.setErrorFunction(
                [this](const Vector &xval, Vector &fval, Matrix &jacobian)
                { this->evaluateErrorFunction(xval, fval, jacobian); });

            Vector xval = initialGuess;
            Vector fval;
            Matrix jacobian;
            Vector gradient;
            Scalar gradLen = minGradLen_ + 1;
            Scalar stepSize;
            Scalar error = minError_ + 1;
            Vector step = Vector::Zero(xval.size());
            Scalar stepLen = minStepLen_ + 1;
            bool callbackResult = true;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                gradLen >= minGradLen_ &&
                stepLen >= minStepLen_ &&
                error >= minError_ &&
                callbackResult)
            {
                xval -= step;
                evaluateErrorFunction(xval, fval, jacobian);
                error = fval.squaredNorm() / 2;
                gradient = jacobian.transpose() * fval;
                gradLen = gradient.norm();

                calculateStep(xval, fval, jacobian, gradient, step);

                // update step according to step size
                stepSize = stepSize_(xval, fval, jacobian, gradient, step);
                step *= stepSize;
                stepLen = step.norm();
                // evaluate callback and save its result
                callbackResult = callback_(iterations + 1, xval, fval, jacobian,
                    gradient, step);

                if(verbosity_ > 0)
                {
                    std::stringstream ss;
                    ss << "it=" << std::setfill('0')
                        << std::setw(4) << iterations
                        << std::fixed << std::showpoint << std::setprecision(6)
                        << "    stepsize=" << stepSize
                        << "    steplen=" << stepLen
                        << "    gradlen=" << gradLen;

                    if(verbosity_ > 1)
                        ss << "    callback=" << (callbackResult ? "true" : "false");

                    ss << "    error=" << error;

                    if(verbosity_ > 2)
                        ss << "    fval=" << vector2str(fval);
                    if(verbosity_ > 3)
                        ss << "    xval=" << vector2str(xval);
                    if(verbosity_ > 4)
                        ss << "    step=" << vector2str(step);
                    std::cout << ss.str() << std::endl;
                }

                ++iterations;
            }

            Result result;
            result.xval = xval;
            result.fval = fval;
            result.error = error;
            result.iterations = iterations;
            result.converged = stepLen < minStepLen_ ||
                gradLen < minGradLen_ ||
                error < minError_;

            return result;
        }

        virtual void calculateStep(const Vector &xval,
            const Vector &fval,
            const Matrix &jacobian,
            const Vector &gradient,
            Vector &step) = 0;

    };

    template<typename Scalar,
        typename ErrorFunction,
        typename StepSize=BarzilaiBorwein<Scalar>,
        typename Callback=NoCallback<Scalar>,
        typename FiniteDifferences=CentralDifferences<Scalar>>
    class GradientDescent : public LeastSquaresAlgorithm<Scalar, ErrorFunction,
        StepSize, Callback, FiniteDifferences>
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;

        GradientDescent()
            : LeastSquaresAlgorithm<Scalar, ErrorFunction,
                StepSize, Callback, FiniteDifferences>()
        { }

        void calculateStep(const Vector &,
            const Vector &,
            const Matrix &,
            const Vector &gradient,
            Vector &step) override
        {
            step = gradient;
        }

    };

    template<typename Scalar,
        typename ErrorFunction,
        typename StepSize=ArmijoBacktracking<Scalar>,
        typename Callback=NoCallback<Scalar>,
        typename FiniteDifferences=CentralDifferences<Scalar>,
        typename Solver=DenseSVDSolver<Scalar>>
    class GaussNewton : public LeastSquaresAlgorithm<Scalar, ErrorFunction,
        StepSize, Callback, FiniteDifferences>
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;

    public:
        GaussNewton()
            : LeastSquaresAlgorithm<Scalar, ErrorFunction,
                StepSize, Callback, FiniteDifferences>()
        { }

        void calculateStep(const Vector &,
            const Vector &,
            const Matrix &jacobian,
            const Vector &gradient,
            Vector &step) override
        {
            Solver solver;

            Matrix A = jacobian.transpose() * jacobian;
            solver(A, gradient, step);
        }
    };

    template<typename Scalar,
        typename ErrorFunction,
        typename Callback=NoCallback<Scalar>,
        typename FiniteDifferences=CentralDifferences<Scalar>,
        typename Solver=DenseSVDSolver<Scalar>>
    class LevenbergMarquardt : public LeastSquaresAlgorithm<Scalar, ErrorFunction,
        ConstantStepSize<Scalar>, Callback, FiniteDifferences>
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
    private:
        Scalar increase_;
        Scalar decrease_;
        Scalar lambda_;
        Index maxItLM_;

    public:
        LevenbergMarquardt()
            : LeastSquaresAlgorithm<Scalar, ErrorFunction,
                ConstantStepSize<Scalar>, Callback, FiniteDifferences>(),
                increase_(static_cast<Scalar>(2)),
                decrease_(static_cast<Scalar>(0.5)),
                lambda_(static_cast<Scalar>(1)),
                maxItLM_(0)
        { }

        /** Set the initial gradient descent factor of levenberg marquardt.
          * @param lambda gradient descent factor */
        void setLambda(const Scalar lambda)
        {
            lambda_ = lambda;
        }

        /** Set maximum iterations of the levenberg marquardt optimization.
          * Set to 0 or negative for infinite iterations.
          * @param iterations maximum iterations for lambda search */
        void setMaxIterationsLM(const Index iterations)
        {
            maxItLM_ = iterations;
        }

        /** Set the increase factor for the lambda damping.
          * Make sure the value is greater than 1.
          * @param increase factor for increasing lambda */
        void setLambdaIncrease(const Scalar increase)
        {
            assert(increase > static_cast<Scalar>(1));
            increase_ = increase;
        }

        /** Set the decrease factor for the lambda damping.
          * Make sure the value is in (0, 1).
          * @param increase factor for increasing lambda */
        void setLambdaDecrease(const Scalar decrease)
        {
            assert(decrease < static_cast<Scalar>(1));
            assert(decrease > static_cast<Scalar>(0));
            decrease_ = decrease;
        }

        void calculateStep(const Vector &xval,
            const Vector &fval,
            const Matrix &jacobian,
            const Vector &gradient,
            Vector &step) override
        {
            Solver solver;
            Scalar error = fval.squaredNorm() / 2;
            Scalar errorN = error + 1;

            Vector xvalN;
            Vector fvalN;
            Matrix jacobianN;

            Matrix jacobianSq = jacobian.transpose() * jacobian;
            Matrix A;

            Index iterations = 0;
            while((maxItLM_ <= 0 || iterations < maxItLM_) &&
                errorN > error)
            {
                A = jacobianSq;
                // add identity matrix
                for(Index i = 0; i < A.rows(); ++i)
                    A(i, i) += lambda_;

                solver(A, gradient, step);

                xvalN = xval - step;
                this->errorFunction_(xvalN, fvalN, jacobianN);
                errorN = fvalN.squaredNorm() / 2;

                if(errorN > error)
                    lambda_ *= increase_;
                else
                    lambda_ *= decrease_;

                ++iterations;
            }
        }
    };

    /** Implementation of Powell's Dogleg Method. */
    template<typename Scalar,
        typename ErrorFunction,
        typename Callback=NoCallback<Scalar>,
        typename FiniteDifferences=CentralDifferences<Scalar>,
        typename Solver=DenseSVDSolver<Scalar>>
    class DoglegMethod : public LeastSquaresAlgorithm<Scalar, ErrorFunction,
        ConstantStepSize<Scalar>, Callback, FiniteDifferences>
    {
    public:
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
        typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;

    private:
        Scalar radius_;
        Scalar maxRadius_;
        Scalar radiusEps_;
        Scalar acceptFitness_;
        Index maxItTR_;

        Scalar calulateModelFitness(const Vector &xval,
            const Vector &fval,
            const Vector &gradient,
            const Matrix &hessian,
            const Vector &step)
        {
            Scalar error = fval.squaredNorm() / 2;

            // evaluate the error function at the new position
            Vector xvalNext = xval + step;
            Vector fvalNext;
            Matrix jacobianNext;
            this->errorFunction_(xvalNext, fvalNext, jacobianNext);
            // compute the actual new error
            Scalar nextError = fvalNext.squaredNorm() / 2;
            // compute the new error by the model
            Scalar modelError = error + gradient.dot(step) + step.dot(hessian * step) / 2;

            return (error - nextError) / (error - modelError);
        }

    public:
        DoglegMethod()
            : LeastSquaresAlgorithm<Scalar, ErrorFunction,
                ConstantStepSize<Scalar>, Callback, FiniteDifferences>(),
                radius_(1),
                maxRadius_(static_cast<Scalar>(2)),
                radiusEps_(static_cast<Scalar>(1e-6)),
                acceptFitness_(static_cast<Scalar>(0.25)),
                maxItTR_(0)
        { }

        /** Set maximum iterations of the trust region radius search.
          * Set to 0 or negative for infinite iterations.
          * @param iterations maximum iterations for radius search */
        void setMaxIterationsTR(const Index iterations)
        {
            maxItTR_ = iterations;
        }

        /** Set the minimum fitness value at which a model is accepted.
          * The model fitness is computed as follows:
          *
          * fitness = f(xval) - f(xval + step) / m(0) - m(step)
          *
          * Where f(x) is the objective error function and m(x) is the
          * model function describe by the trust region method.
          *
          * @param fitness minimum fitness for step acceptance */
        void setAcceptanceFitness(const Scalar fitness)
        {
            acceptFitness_ = fitness;
        }

        /** Set the comparison epsilon on how close the step should be
          * to the trust region radius to trigger an increase of the radius.
          * Should usually be picked low, e.g. 1e-8.
          * @param eps comparison epsilon for radius increase */
        void setRaidusEps(const Scalar eps)
        {
            radiusEps_ = eps;
        }

        void calculateStep(const Vector &xval,
            const Vector &fval,
            const Matrix &jacobian,
            const Vector &gradient,
            Vector &step) override
        {
            // approximate hessian
            Matrix hessian = jacobian.transpose() * jacobian;

            // compute the full newton step
            Vector fullStep;
            Solver solver;
            solver(hessian, gradient, fullStep);
            fullStep = -fullStep;

            // precompute the full step length
            Scalar fullStepLenSq = fullStep.squaredNorm();

            // compute the cauchy step
            Scalar gradientLenSq = gradient.squaredNorm();
            Scalar curvature = gradient.dot(hessian * gradient);
            Vector cauchyStep = -(gradientLenSq / curvature) * gradient;
            Scalar cauchyStepLenSq = cauchyStep.squaredNorm();

            // compute step diff
            Vector diffStep = fullStep - cauchyStep;
            Scalar diffLenSq = diffStep.squaredNorm();
            Scalar diffFac = cauchyStep.dot(diffStep) / diffLenSq;

            Scalar modelFitness = acceptFitness_ - 1;
            Index iteration = 0;

            // keep computing while the model fitness is bad
            while(modelFitness < acceptFitness_
                && (maxItTR_ <= 0 || iteration < maxItTR_))
            {
                Scalar radiusSq = radius_ * radius_;

                // if the full step is within the trust region simply
                // use it, it provides a good minimizer
                if(fullStepLenSq <= radiusSq)
                {
                    step = fullStep;
                }
                else
                {
                    // if the cauchy step lies outside the trust region
                    // go towards it until the trust region boundary
                    if(cauchyStepLenSq >= radiusSq)
                    {
                        step = (radius_ / std::sqrt(cauchyStepLenSq)) * cauchyStep;
                    }
                    else
                    {
                        Scalar secondTerm = std::sqrt(diffFac * diffFac + (radiusSq + cauchyStepLenSq) / diffLenSq);
                        Scalar scale1 = -diffFac - secondTerm;
                        Scalar scale2 = -diffFac + secondTerm;

                        step = cauchyStep + std::max(scale1, scale2) * (fullStep - cauchyStep);
                    }
                }

                // compute the model fitness to determine the update scheme for
                // the trust region radius
                modelFitness = calulateModelFitness(xval, fval, gradient, hessian, step);

                Scalar stepLen = step.norm();

                // if the model fitness is really bad reduce the radius!
                if(modelFitness < static_cast<Scalar>(0.25))
                    radius_ = static_cast<Scalar>(0.25) * stepLen;
                // if the model fitness is very good then increase it
                else if(modelFitness > static_cast<Scalar>(0.75) && std::abs(stepLen - radius_) < radiusEps_)
                {
                    // use the double radius
                    radius_ = 2 * radius_;
                    // maintain radius border if configured
                    if(maxRadius_ > 0 && radius_ > maxRadius_)
                        radius_ = maxRadius_;
                }

                ++iteration;
            }

            step = -step;
        }
    };
}

#endif
