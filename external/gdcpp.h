/* gdcpp.h
 *
 *     Author: Fabian Meyer
 * Created On: 12 Jul 2019
 *    License: MIT
 */

#ifndef GDCPP_GDCPP_H_
#define GDCPP_GDCPP_H_

#include <Eigen/Geometry>
#include <limits>
#include <iostream>
#include <iomanip>
#include <functional>

namespace gdc
{
    typedef long int Index;

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
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &)> Objective;
    private:
        Scalar eps_;
        Index threads_;
        Objective objective_;
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

        void setThreads(const Index threads)
        {
            threads_ = threads;
        }

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Scalar fval,
            Vector &gradient)
        {
            assert(objective_);

            gradient.resize(xval.size());
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < xval.size(); ++i)
            {
                Vector xvalN = xval;
                xvalN(i) += eps_;
                Scalar fvalN = objective_(xvalN);

                gradient(i) = (fvalN - fval) / eps_;
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
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &)> Objective;
    private:
        Scalar eps_;
        Index threads_;
        Objective objective_;
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

        void setThreads(const Index threads)
        {
            threads_ = threads;
        }

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Scalar fval,
            Vector &gradient)
        {
            assert(objective_);

            gradient.resize(xval.size());
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < xval.size(); ++i)
            {
                Vector xvalN = xval;
                xvalN(i) -= eps_;
                Scalar fvalN = objective_(xvalN);
                gradient(i) = (fval - fvalN) / eps_;
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
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &)> Objective;
    private:
        Scalar eps_;
        Index threads_;
        Objective objective_;
    public:
        CentralDifferences()
            : CentralDifferences(
                gismo::math::sqrt(std::numeric_limits<Scalar>::epsilon()))
        { }

        CentralDifferences(const Scalar eps)
            : eps_(eps), threads_(1), objective_()
        { }

        void setNumericalEpsilon(const Scalar eps)
        {
            eps_ = eps;
        }

        void setThreads(const Index threads)
        {
            threads_ = threads;
        }

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void operator()(const Vector &xval,
            const Scalar,
            Vector &gradient)
        {
            assert(objective_);

            Vector fvals(xval.size() * 2);
            #pragma omp parallel for num_threads(threads_)
            for(Index i = 0; i < fvals.size(); ++i)
            {
                Index idx = i / 2;
                Vector xvalN = xval;
                if(i % 2 == 0)
                    xvalN(idx) += eps_ / 2;
                else
                    xvalN(idx) -= eps_ / 2;

                fvals(i) = objective_(xvalN);
            }

            gradient.resize(xval.size());
            for(Index i = 0; i < xval.size(); ++i)
                gradient(i) = (fvals(i * 2) - fvals(i * 2 + 1)) / eps_;
        }
    };

    template<typename Scalar>
    struct NoDerivativeApproximation
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &)> Objective;
    public:
        NoDerivativeApproximation()
        { }

        NoDerivativeApproximation(const Scalar)
        { }

        void setNumericalEpsilon(const Scalar)
        { }

        void setThreads(const Index)
        { }

        void setObjective(const Objective &)
        { }

        void operator()(const Vector &,
                        const Scalar,
                        Vector &)
        { }
    };

    /** Dummy callback functor, which does nothing. */
    template<typename Scalar>
    struct NoCallback
    {
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

        bool operator()(const Index,
            const Vector &,
            const Scalar,
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
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &, Vector &)> Objective;
        typedef std::function<void(const Vector &, const Scalar, Vector &)> FiniteDifferences;
    private:
        Scalar stepSize_;
    public:

        ConstantStepSize()
            : ConstantStepSize(0.7)
        { }

        ConstantStepSize(const Scalar stepSize)
            : stepSize_(stepSize)
        { }

        /** Set the step size returned by this functor.
          * @param stepSize step size returned by functor */
        void setStepSize(const Scalar stepSize)
        {
            stepSize_ = stepSize;
        }

        void setObjective(const Objective &)
        { }

        void setFiniteDifferences(const FiniteDifferences &)
        { }

        Scalar operator()(const Vector &,
            const Scalar,
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
      * y_k = grad_k - grad_k-1   k >= 1
      * Direct:  stepSize = (s_k^T * s_k) / (y_k^T * s_k)
      * Inverse: stepSize = (y_k^T * s_k) / (y_k^T * y_k)
      *
      * The very first step is computed as a constant. */
    template<typename Scalar>
    class BarzilaiBorwein
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &, Vector &)> Objective;
        typedef std::function<void(const Vector &, const Scalar, Vector &)> FiniteDifferences;

        enum class Method
        {
            Direct,
            Inverse
        };
    private:
        Vector lastXval_;
        Vector lastGradient_;
        Method method_;
        Scalar constStep_;

        Scalar constantStep() const
        {
            return constStep_;
        }

        Scalar directStep(const Vector &xval,
            const Vector &gradient)
        {
            Vector sk = xval - lastXval_;
            Vector yk = gradient - lastGradient_;
            Scalar num = sk.dot(sk);
            Scalar denom = sk.dot(yk);

            if(denom == 0)
                return 1;
            else
                return gismo::math::abs(num / denom);
        }

        Scalar inverseStep(const Vector &xval,
            const Vector &gradient)
        {
            Vector sk = xval - lastXval_;
            Vector yk = gradient - lastGradient_;
            Scalar num = sk.dot(yk);
            Scalar denom = yk.dot(yk);

            if(denom == 0)
                return 1;
            else
                return gismo::math::abs(num / denom);
        }
    public:
        BarzilaiBorwein()
            : BarzilaiBorwein(Method::Direct, 1e-4)
        { }

        BarzilaiBorwein(const Method method, const Scalar constStep)
            : lastXval_(), lastGradient_(), method_(method),
            constStep_(constStep)
        { }

        void setObjective(const Objective &)
        { }

        void setFiniteDifferences(const FiniteDifferences &)
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
            const Scalar,
            const Vector &gradient)
        {
            Scalar stepSize = 0;
            if(lastXval_.size() == 0)
            {
                stepSize = constStep_;
            }
            else
            {
                switch(method_)
                {
                case Method::Direct:
                    stepSize = directStep(xval, gradient);
                    break;
                case Method::Inverse:
                    stepSize = inverseStep(xval, gradient);
                    break;
                default:
                    assert(false);
                    break;
                }
            }

            lastGradient_ = gradient;
            lastXval_ = xval;

            return stepSize;
        }
    };


    /** Step size functor to perform Armijo Linesearch with backtracking.
      * The functor iteratively decreases the step size until the following
      * conditions are met:
      *
      * Armijo: f(x - stepSize * grad(x)) <= f(x) - cArmijo * stepSize * grad(x)^T * grad(x)
      *
      * If either condition does not hold the step size is decreased:
      *
      * stepSize = decrease * stepSize */
    template<typename Scalar>
    class ArmijoBacktracking
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &, Vector &)> Objective;
        typedef std::function<void(const Vector &, const Scalar, Vector &)> FiniteDifferences;
    protected:
        Scalar decrease_;
        Scalar cArmijo_;
        Scalar minStep_;
        Scalar maxStep_;
        Index maxIt_;
        Objective objective_;
        FiniteDifferences finiteDifferences_;

        Scalar evaluateObjective(const Vector &xval, Vector &gradient)
        {
            gradient.resize(0);
            Scalar fval = objective_(xval, gradient);
            if(gradient.size() == 0)
                finiteDifferences_(xval, fval, gradient);
            return fval;
        }

        virtual bool computeSecondCondition(const Scalar,
            const Scalar,
            const Scalar,
            const Vector &,
            const Vector &)
        {
            return true;
        }
    public:
        ArmijoBacktracking()
            : ArmijoBacktracking(0.8, 1e-4, 1e-12, 1.0, 0)
        { }

        ArmijoBacktracking(const Scalar decrease,
            const Scalar cArmijo,
            const Scalar minStep,
            const Scalar maxStep,
            const Index iterations)
            : decrease_(decrease), cArmijo_(cArmijo), minStep_(minStep),
            maxStep_(maxStep), maxIt_(iterations), objective_()
        {
            assert(decrease > 0);
            assert(decrease < 1);
            assert(cArmijo > 0);
            assert(cArmijo < 0.5);
            assert(minStep < maxStep);
        }

        /** Set the decreasing factor for backtracking.
          * Assure that decrease in (0, 1).
          * @param decrease decreasing factor */
        void setBacktrackingDecrease(const Scalar decrease)
        {
            assert(decrease > 0);
            assert(decrease < 1);
            decrease_ = decrease;
        }

        /** Set the relaxation constant for the Armijo condition (see class
          * description).
          * Assure cArmijo in (0, 0.5).
          * @param cArmijo armijo constant */
        void setArmijoConstant(const Scalar cArmijo)
        {
            assert(cArmijo > 0);
            assert(cArmijo < 0.5);
            cArmijo_ = cArmijo;
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

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void setFiniteDifferences(const FiniteDifferences &finiteDifferences)
        {
            finiteDifferences_ = finiteDifferences;
        }

        Scalar operator()(const Vector &xval,
            const Scalar fval,
            const Vector &gradient)
        {
            assert(objective_);
            assert(finiteDifferences_);

            Scalar stepSize = maxStep_ / decrease_;
            Vector gradientN;
            Vector xvalN;
            Scalar fvalN;
            bool armijoCondition = false;
            bool secondCondition = false;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                stepSize * decrease_ >= minStep_ &&
                !(armijoCondition && secondCondition))
            {
                stepSize = decrease_ * stepSize;
                xvalN = xval - stepSize * gradient;
                fvalN = evaluateObjective(xvalN, gradientN);

                armijoCondition = fvalN <= fval - cArmijo_ * stepSize * gradient.dot(gradient);
                secondCondition = computeSecondCondition(stepSize, fval, fvalN, gradient, gradientN);

                ++iterations;
            }

            return stepSize;
        }
    };

    /** Step size functor to perform Wolfe Linesearch with backtracking.
      * The functor iteratively decreases the step size until the following
      * conditions are met:
      *
      * Armijo: f(x - stepSize * grad(x)) <= f(x) - cArmijo * stepSize * grad(x)^T * grad(x)
      * Wolfe: grad(x)^T grad(x - stepSize * grad(x)) <= cWolfe * grad(x)^T * grad(x)
      *
      * If either condition does not hold the step size is decreased:
      *
      * stepSize = decrease * stepSize */
    template<typename Scalar>
    class WolfeBacktracking : public ArmijoBacktracking<Scalar>
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &, Vector &)> Objective;
        typedef std::function<void(const Vector &, const Scalar, Vector &)> FiniteDifferences;
    protected:
        Scalar cWolfe_;

        virtual bool computeSecondCondition(const Scalar,
            const Scalar ,
            const Scalar ,
            const Vector &gradient,
            const Vector &gradientN)
        {
            return gradient.dot(gradientN) <= cWolfe_ * gradient.dot(gradient);
        }
    public:
        WolfeBacktracking()
            : WolfeBacktracking(0.8, 1e-4, 0.9, 1e-12, 1.0, 0)
        { }

        WolfeBacktracking(const Scalar decrease,
            const Scalar cArmijo,
            const Scalar cWolfe,
            const Scalar minStep,
            const Scalar maxStep,
            const Index iterations)
            : ArmijoBacktracking<Scalar>(decrease, cArmijo, minStep, maxStep,
                iterations),cWolfe_(cWolfe)
        {
            assert(cWolfe < 1);
            assert(cArmijo < cWolfe);
        }

        /** Set the wolfe constants for Armijo and Wolfe condition (see class
          * description).
          * Assure that c1 < c2 < 1 and c1 in (0, 0.5).
          * @param c1 armijo constant
          * @param c2 wolfe constant */
        void setWolfeConstant(const Scalar cWolfe)
        {
            assert(cWolfe < 1);
            cWolfe_ = cWolfe;
        }
    };

    /** Step size functor which searches for a step that reduces the function
      * value.
      * The functor iteratively decreases the step size until the following
      * condition is met:
      *
      * f(x - stepSize * grad) < f(x)
      *
      * If this condition does not hold the step size is decreased:
      *
      * stepSize = decrease * stepSize
      *
      * This functor does not require to compute any gradients and does not use
      * finite differences. */
    template<typename Scalar>
    class DecreaseBacktracking
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef std::function<Scalar(const Vector &, Vector &)> Objective;
        typedef std::function<void(const Vector &, const Scalar, Vector &)> FiniteDifferences;
    private:
        Scalar decrease_;
        Scalar minStep_;
        Scalar maxStep_;
        Index maxIt_;
        Objective objective_;

    public:
        DecreaseBacktracking()
            : DecreaseBacktracking(0.8, 1e-12, 1.0, 0)
        { }

        DecreaseBacktracking(const Scalar decrease,
            const Scalar minStep,
            const Scalar maxStep,
            const Index iterations)
            : decrease_(decrease), minStep_(minStep),
            maxStep_(maxStep), maxIt_(iterations), objective_()
        { }

        /** Set the decreasing factor for backtracking.
          * Assure that decrease in (0, 1).
          * @param decrease decreasing factor */
        void setBacktrackingDecrease(const Scalar decrease)
        {
            decrease_ = decrease;
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

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void setFiniteDifferences(const FiniteDifferences &)
        { }

        Scalar operator()(const Vector &xval,
            const Scalar fval,
            const Vector &gradient)
        {
            assert(objective_);

            Scalar stepSize = maxStep_ / decrease_;
            Vector xvalN;
            Vector gradientN;
            Scalar fvalN;
            bool improvement = false;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                stepSize * decrease_ >= minStep_ &&
                !improvement)
            {
                stepSize = decrease_ * stepSize;
                xvalN = xval - stepSize * gradient;
                fvalN = objective_(xvalN, gradientN);

                improvement = fvalN < fval;

                ++iterations;
            }

            return stepSize;
        }
    };

    template<typename Scalar,
        typename Objective,
        typename StepSize=BarzilaiBorwein<Scalar>,
        typename Callback=NoCallback<Scalar>,
        typename FiniteDifferences=CentralDifferences<Scalar> >
    class GradientDescent
    {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

        struct Result
        {
            Index iterations;
            bool converged;
            Scalar fval;
            Vector xval;
        };
    private:
        /*
         template <typename U>
         static constexpr auto hasGradient(U * t) ->
             decltype(t->operator Scalar(std::declval<const Vector>(), std::declval<Vector>()), std::true_type{});

        template <typename U>
        static constexpr auto hasGradient(...) -> std::false_type;
        */
    protected:
        Index maxIt_;
        Scalar minGradientLen_;
        Scalar minStepLen_;
        Scalar momentum_;
        Index verbosity_;
        Objective objective_;
        StepSize stepSize_;
        Callback callback_;
        FiniteDifferences finiteDifferences_;
        
//        typename std::enable_if<decltype(hasGradient<Objective>(nullptr)),Scalar>::type
        Scalar evaluateObjective(const Vector &xval, Vector &gradient)
        {
            gradient.resize(0);
            Scalar fval = objective_(xval, gradient);
            if(gradient.size() == 0)
                finiteDifferences_(xval, fval, gradient);
            return fval;
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

        GradientDescent()
            : maxIt_(0), minGradientLen_(static_cast<Scalar>(1e-9)),
            minStepLen_(static_cast<Scalar>(1e-9)), momentum_(0),
            verbosity_(0), objective_(), stepSize_(), callback_(),
            finiteDifferences_()
        {

        }

        ~GradientDescent()
        {

        }

        void setThreads(const Index threads)
        {
            finiteDifferences_.setThreads(threads);
        }

        void setNumericalEpsilon(const Scalar eps)
        {
            finiteDifferences_.setNumericalEpsilon(eps);
        }

        void setMaxIterations(const Index iterations)
        {
            maxIt_ = iterations;
        }

        void setObjective(const Objective &objective)
        {
            objective_ = objective;
        }

        void setCallback(const Callback &callback)
        {
            callback_ = callback;
        }

        void setMinGradientLength(const Scalar gradientLen)
        {
            minGradientLen_ = gradientLen;
        }

        void setMinStepLength(const Scalar stepLen)
        {
            minStepLen_ = stepLen;
        }

        void setStepSize(const StepSize stepSize)
        {
            stepSize_ = stepSize;
        }

        void setMomentum(const Scalar momentum)
        {
            momentum_ = momentum;
        }

        void setVerbosity(const Index verbosity)
        {
            verbosity_ = verbosity;
        }

        Result minimize(const Vector &initialGuess)
        {
            finiteDifferences_.setObjective(
                [this](const Vector &xval)
                { Vector tmp; return this->objective_(xval, tmp); });
            stepSize_.setObjective(
                [this](const Vector &xval, Vector &gradient)
                { return this->objective_(xval, gradient); });
            stepSize_.setFiniteDifferences(
                [this](const Vector &xval, const Scalar fval, Vector &gradient)
                { this->finiteDifferences_(xval, fval, gradient); });

            Vector xval = initialGuess;
            Vector gradient;
            Scalar fval(0);
            Scalar gradientLen = minGradientLen_ + 1;
            Scalar stepSize;
            Vector step = Vector::Zero(xval.size());
            Scalar stepLen = minStepLen_ + 1;
            bool callbackResult = true;

            Index iterations = 0;
            while((maxIt_ <= 0 || iterations < maxIt_) &&
                gradientLen >= minGradientLen_ &&
                stepLen >= minStepLen_
                && callbackResult)
            {
                xval -= step;
                fval = evaluateObjective(xval, gradient);
                gradientLen = gradient.norm();
                // update step according to step size and momentum
                stepSize = stepSize_(xval, fval, gradient);
                step = momentum_ * step + (1 - momentum_) * stepSize * gradient;
                stepLen = step.norm();
                // evaluate callback an save its result
                callbackResult = callback_(iterations, xval, fval, gradient);

                if(verbosity_ > 0)
                {
                    std::stringstream ss;
                    ss << "it=" << std::setfill('0')
                        << std::setw(4) << iterations
                        << std::fixed << std::showpoint << std::setprecision(6)
                        << "    gradlen=" << gradientLen
                        << "    stepsize=" << stepSize
                        << "    steplen=" << stepLen;

                    if(verbosity_ > 2)
                        ss << "    callback=" << (callbackResult ? "true" : "false");

                    ss << "    fval=" << fval;

                    if(verbosity_ > 1)
                        ss << "    xval=" << vector2str(xval);
                    if(verbosity_ > 2)
                        ss << "    gradient=" << vector2str(gradient);
                    if(verbosity_ > 3)
                        ss << "    step=" << vector2str(step);
                    std::cout << ss.str() << std::endl;
                }

                ++iterations;
            }

            Result result;
            result.xval = xval;
            result.fval = fval;
            result.iterations = iterations;
            result.converged = gradientLen < minGradientLen_ ||
                stepLen < minStepLen_;

            return result;
        }
    };
}

#endif
