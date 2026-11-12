# SUPPLEMENTARY TERMS FOR ARTIFICIAL INTELLIGENCE USE

*Additional Terms under Section 7 of the GNU AGPL version 3*  
**Version 1.3 — Effective November 13, 2026**

These supplementary terms apply to the Plasmid source code and associated materials (the "Program"). They are issued pursuant to Section 7 of the GNU Affero General Public License, Version 3 ("AGPL-3.0"), under which the Program is licensed. 

These terms are additional restrictions permitted under Section 7 of the AGPL-3.0. They supplement, and do not replace, the AGPL-3.0. In the event of a conflict between these terms and the AGPL-3.0, the AGPL-3.0 shall govern, except to the extent that these terms are valid additional terms under Section 7 of the AGPL-3.0.

---

## 1. Conditional Restriction on Artificial Intelligence Training Data

Except as explicitly provided in Section 2 of these supplementary terms, the Program, including its source code, documentation, and associated materials (collectively, the "Materials"), may not be used, in whole or in part, as training data, fine-tuning data, reinforcement learning data, evaluation data, or any other form of input for machine learning models, neural networks, large language models, or other artificial intelligence systems, without the prior express written permission of the copyright holder of Plasmid.

This restriction applies irrespective of the method of ingestion, including but not limited to:

- Direct copying of source files into training pipelines;
- Automated scraping or crawling of the repository;
- Inclusion in curated datasets that aggregate open-source code;
- Embedding or vectorizing source code for retrieval-augmented generation (RAG) systems.

> **Clarification:** This restriction does not prohibit individual developers from utilizing general-purpose AI coding assistants (such as GitHub Copilot, Cursor, or analogous tools) in their ordinary development workflow when contributing to this Program. The restriction is expressly directed at the systematic, large-scale use of this Program's source code as input for training or improving artificial intelligence systems.

## 2. Copyleft Extension to AI Model Artifacts

In the event that the Materials are utilized, in whole or in part, as training data, fine-tuning data, or any other form of input for machine learning models, neural networks, large language models, or other artificial intelligence systems (whether pursuant to written permission or otherwise), the resulting artificial intelligence model shall be deemed a "Modified Version" or "Derived Work" of the Program under the AGPL-3.0.

Accordingly, the model developer or any entity conveying such a model is strictly obligated to comply with all provisions of the AGPL-3.0, including the requirement to provide the complete "Corresponding Source" to all recipients. For the avoidance of doubt, the "Corresponding Source" for such an AI-derived work explicitly encompasses, but is not limited to:

- The complete source code of the resulting model, including all inference, deployment, and orchestration scripts;
- The trained model weights and parameters;
- The complete training methodology, including data pipelines, preprocessing steps, hyperparameters, and training scripts;
- Any prompt engineering strategies, agent orchestration patterns, or design elements derived from the Materials.

Any failure to make these artifacts publicly available under the terms of the AGPL-3.0 constitutes a material breach of these supplementary terms and the underlying license, resulting in the automatic termination of your rights under the AGPL-3.0 as specified in Section 8 of the License.

## 3. Permitted AI-Assisted Reproduction and Reverse Engineering

The use of artificial intelligence tools, including but not limited to large language models, code generation systems, or automated analysis tools, to analyze, reverse-engineer, or reproduce the architecture, algorithms, data structures, prompt engineering strategies, agent orchestration patterns, or design patterns of this Program is expressly permitted. 

However, any such reproduction or reverse engineering is strictly conditional. If you utilize artificial intelligence tools to create a work based on, derived from, or substantially replicating the Program (a "Derived AI Work"), you must comply with the following obligations:

1. **Open Source Requirement:** The complete source code of the Derived AI Work must be made publicly available and distributed as open-source software.
2. **AGPL-3.0 Licensing:** The Derived AI Work must be licensed to all downstream recipients under the terms of the GNU Affero General Public License, Version 3 (AGPL-3.0).
3. **Retention of Supplementary Terms:** The Derived AI Work must be accompanied by, and explicitly subject to, these *Supplementary Terms for Artificial Intelligence Use* in their entirety. All restrictions and obligations herein, including but not limited to the AI training restrictions (Section 1) and copyleft extensions (Section 2), shall apply with equal force to the Derived AI Work.

> **Clarification:** This section ensures that AI-assisted development remains a collaborative, open-source endeavor. It prevents the proprietary enclosure of the Program's core design while fully permitting innovation, adaptation, and the creation of competing or similar products, provided they remain within the AGPL-3.0 copyleft ecosystem.

## 4. Legal Basis Under AGPL-3.0 Section 7

These supplementary terms are issued under the authority of Section 7 of the AGPL-3.0, specifically:

- **Section 7(b)**, which permits requiring the preservation of specified reasonable legal notices or author attributions in the Materials and in any works derived from the Materials, including through AI-mediated derivation;
- **Section 7(c)**, which prohibits misrepresentation of the origin of the Materials. The use of artificial intelligence tools to obscure the origin of code, weights, or designs derived from this Program constitutes misrepresentation of origin within the meaning of this section.

## 5. Severability

If any provision of these supplementary terms is held to be unenforceable or invalid by a court of competent jurisdiction, that provision shall be modified to the minimum extent necessary to make it enforceable, or, if modification is not possible, severed from these terms. The remaining provisions shall continue in full force and effect.

---

## Appendix: How to Apply These Supplementary Terms to Your Files

To apply these supplementary terms to your work, attach the following boilerplate notice to the start of each source file. The "copyright" line and a pointer to where the full notice is found should be included:

```text
Copyright (C) [2026] [Zhewei Chen]

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version, subject to the Supplementary Terms 
for Artificial Intelligence Use (Version 1.3, Effective November 13, 2026).

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
