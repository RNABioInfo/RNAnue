#pragma once

// Internal
#include "AnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

class EvaluationMetrics {
   public:
    EvaluationMetrics(double pvalue, double padj = std::numeric_limits<double>::quiet_NaN())
        : pvalue(pvalue), padj(padj) {}

    [[nodiscard]] auto getPValue() const -> double { return pvalue; }
    [[nodiscard]] auto getPadj() const -> double { return padj; }

    void setPadj(double padj) { this->padj = padj; }

   private:
    double pvalue;
    double padj;
};

class EvaluatedInteractionCluster : public AnnotatedInteractionCluster {
   public:
    EvaluatedInteractionCluster(AnnotatedInteractionCluster &&cluster, double pvalue,
                                double padj = std::numeric_limits<double>::quiet_NaN())
        : AnnotatedInteractionCluster(std::move(cluster)), metrics(pvalue, padj) {}

    EvaluatedInteractionCluster(const AnnotatedInteractionCluster &cluster, double pvalue,
                                double padj = std::numeric_limits<double>::quiet_NaN())
        : AnnotatedInteractionCluster(cluster), metrics(pvalue, padj) {}

    [[nodiscard]] auto getPValue() const -> double { return metrics.getPValue(); }
    [[nodiscard]] auto getPadj() const -> double { return metrics.getPadj(); }

    void setPadj(double padj) { metrics.setPadj(padj); }

   private:
    EvaluationMetrics metrics;
};

}  // namespace pipelines::analyze
